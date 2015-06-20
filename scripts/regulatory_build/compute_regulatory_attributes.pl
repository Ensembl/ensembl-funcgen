#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=head1 NAME

load_build.pl

=head1 SYNOPSIS

perl load_new_assembly.pl --host $host --user $user --pass $pass --dbname homo_sapiens_funcgen_76_38 --base_dir hg38/

=head1 DESCRIPTION

Loads a segmentation BigBed file annotated by the new regulatory build into the database.
Params:
	* base_dir: directory with assembly name (e.g. ./hg38) create by build
	* your usual Ensembl Funcgen MySQL params: host, user etc...

In short, the file $base_dir/segmentations/$segmentation/$cell_type_name.bb must exist

Also bigBedToBed must be on the commandline.

This script is pretty database intensive (think inserting 4Mo entries) so keep an eye 
when running in parallel. In the good times, each job takes ~1h to run. 

=cut

use strict;
use Getopt::Long;
use File::Basename;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::Helper;
use File::Temp qw/ tempfile tempdir /;
${File::Temp::KEEP_ALL} = 1;
our $start_time = time;

main();

####################################################
## Overview
####################################################
sub main {
  print_log("Getting options\n");
  my $options = get_options();
  print_log("Creating regulatory_annotation table\n");
  compute_regulatory_annotations($options);
}

########################################################
## System calls 
## Wrapper function for system calls
## Params:
## - Command line
## Actions:
## - Runs command, prints out error in case of failure
########################################################

sub run {
  my ($cmd) = @_;
  print_log("Running $cmd\n");
  my $exit_code = system($cmd);
  if ($exit_code != 0) {
    die("Failure when running command\n$cmd\n")
  }
}

########################################################
## Print conveninence
## Params:
## - String
## Actions:
## - Prints string, with time stamp in front
########################################################

sub print_log {
  my ($str) = @_;
  my $runtime = time - $start_time;
  print "[$runtime] $str";
}

####################################################
## Command line options
####################################################
sub get_options {
  my %options = ();
  GetOptions (
    \%options,
    "pass|p=s",
    "port=s",
    "host|h=s",
    "user|u=s",
    "dbname|d=s",
  ) or pod2usage( -exitval => 1);
  defined $options{host} || die ("You must define the destination host!\t--host XXXX\n");
  defined $options{user} || die ("You must define the user login!\t--user XXXX\n");
  defined $options{dbname} || die ("You must define the database name!\t--dbname XXXX\n");
  return \%options;
}

#####################################################
# Assign motifs and annotations to regulatory features
#####################################################
# Params: - host
#         - user
#         - pass
#         - dbname
#####################################################

sub compute_regulatory_annotations {
  my ($options) = @_;

  my ($tmp_fh, $cell_type_regulatory_features) = tempfile();
  my ($tmp_fh, $multicell_regulatory_features) = tempfile();
  my ($tmp_fh2, $annotations) = tempfile();
  my ($tmp_fh3, $motifs) = tempfile();
  my ($tmp_fh4, $out) = tempfile();
  close $tmp_fh;
  close $tmp_fh2;
  close $tmp_fh3;
  close $tmp_fh4;

  # Extract regulatory, annotated, and motif features into flat tab-delimited files.
  # The common format of these files is:
  # chrom	start	end	ID
  run("mysql --quick -NB -h $options->{host} -u $options->{user} -p$options->{pass} -D $options->{dbname} -e 'select rf.seq_region_id, rf.seq_region_start - rf.bound_start_length, rf.seq_region_end + rf.bound_end_length, rf.regulatory_feature_id, fs.cell_type_id FROM regulatory_feature rf  JOIN feature_set fs USING(feature_set_id)  JOIN cell_type USING(cell_type_id) WHERE cell_type.name LIKE \"MultiCell\" ORDER by rf.seq_region_id, rf.seq_region_start - rf.bound_start_length'| sort -k1,1 -k2,2n > $multicell_regulatory_features");
  run("mysql --quick -NB -h $options->{host} -u $options->{user} -p$options->{pass} -D $options->{dbname} -e 'select rf.seq_region_id, rf.seq_region_start - rf.bound_start_length, rf.seq_region_end + rf.bound_end_length, rf.regulatory_feature_id, fs.cell_type_id FROM regulatory_feature rf  JOIN feature_set fs USING(feature_set_id)  JOIN cell_type USING(cell_type_id) WHERE cell_type.name NOT LIKE \"MultiCell\" ORDER by rf.seq_region_id, rf.seq_region_start - rf.bound_start_length'| sort -k1,1 -k2,2n > $cell_type_regulatory_features");
  run("mysql --quick -NB -h $options->{host} -u $options->{user} -p$options->{pass} -D $options->{dbname} -e 'select af.seq_region_id, af.seq_region_start, af.seq_region_end, af.annotated_feature_id, fs.cell_type_id FROM annotated_feature af  JOIN feature_set fs USING(feature_set_id) ORDER BY af.seq_region_id, af.seq_region_start;' | awk '\$2 <= \$3' > $annotations");
  run("mysql --quick -NB -h $options->{host} -u $options->{user} -p$options->{pass} -D $options->{dbname} -e 'select mf.seq_region_id, mf.seq_region_start, mf.seq_region_end, mf.motif_feature_id, fs.cell_type_id from motif_feature mf  JOIN associated_motif_feature amf USING(motif_feature_id)  JOIN annotated_feature af USING(annotated_feature_id)  JOIN feature_set fs USING(feature_set_id) GROUP BY motif_feature_id, cell_type_id ORDER BY seq_region_id, seq_region_start;' > $motifs");

  # Overlap regulatory features with (annotated|motif) features. Store ID pairs into one flat file
  run("bedtools intersect -wa -wb -a $cell_type_regulatory_features -b $annotations | awk 'BEGIN {OFS = \"\\t\"} \$5==\$10 {print \$4,\$9,\"annotated\"} ' > $out");
  run("bedtools intersect -wa -wb -a $multicell_regulatory_features -b $motifs | awk 'BEGIN {OFS = \"\\t\"} {print \$4,\$9,\"motif\"} ' >> $out");
  run("bedtools intersect -wa -wb -a $cell_type_regulatory_features -b $motifs | awk 'BEGIN {OFS = \"\\t\"} \$5==\$10 {print \$4,\$9,\"motif\"} ' >> $out");

  # Load into database
  run("mysql -h $options->{host} -u $options->{user} -p$options->{pass} -D $options->{dbname} -e 'TRUNCATE TABLE regulatory_attribute;'");
  run("mysql -h $options->{host} -u $options->{user} -p$options->{pass} -D $options->{dbname} -e 'LOAD DATA LOCAL INFILE \"$out\" INTO TABLE regulatory_attribute;'");

  unlink $cell_type_regulatory_features;
  unlink $multicell_regulatory_features;
  unlink $annotations;
  unlink $motifs;
  unlink $out;
}

1;
