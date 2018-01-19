#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

perl compute_regulatory_attributes.pl --host $host --user $user --pass $pass --dbname homo_sapiens_funcgen_76_38

=head1 DESCRIPTION

This scripts populates the regualaltory_attributes table by dumping regulatory_feature, annotated_feature 
and motif_feature records to bed files. Bedtools is then used to perform intersects which to identify 
overlapping features, including those which overhang the ends of a given regulatory feature.

=cut

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd );
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
  print_log("Creating regulatory_attributes table\n");
  compute_regulatory_attributes($options);
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
  my %options;
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

sub compute_regulatory_attributes {
  my $options = shift;

  my ($tmp_fh,  $epigenome_regulatory_features) = tempfile();
#   my ($tmp_fh1, $multicell_regulatory_features) = tempfile();
  my ($tmp_fh2, $annotations) = tempfile();
  my ($tmp_fh3, $motifs) = tempfile();
  my ($tmp_fh4, $out) = tempfile();
  close $tmp_fh;
#   close $tmp_fh1;
  close $tmp_fh2;
  close $tmp_fh3;
  close $tmp_fh4;

  # Extract regulatory, annotated, and motif features into flat tab-delimited files.
  # The common format of these files is:
  # chrom	start	end	ID 

  run_system_cmd(
    "mysql --quick -NB -h $options->{host} -u $options->{user} -p$options->{pass} -D $options->{dbname} "
    . q( -e "
    select 
      rf.seq_region_id, 
      rf.seq_region_start - rf.bound_start_length as start, 
      rf.seq_region_end + rf.bound_end_length, 
      rf.regulatory_feature_id, 
      ra.epigenome_id 
    FROM 
      regulatory_feature rf
    JOIN 
      regulatory_activity ra 
    USING 
      (regulatory_feature_id)
    JOIN 
      epigenome USING (epigenome_id) 
    ORDER by 
      rf.seq_region_id, start" ) . " > $epigenome_regulatory_features",
     undef, 1); 

  run_system_cmd("mysql --quick -NB -h $options->{host} -u $options->{user} -p$options->{pass} -D $options->{dbname} ".
    '-e "select af.seq_region_id, af.seq_region_start, af.seq_region_end, af.annotated_feature_id, fs.epigenome_id '.
    'FROM annotated_feature af JOIN feature_set fs USING(feature_set_id) where af.seq_region_start <= af.seq_region_end ORDER BY af.seq_region_id, af.seq_region_start" '.
    " > $annotations",  # filter out features where start > end. These should have been caught by now in the SeqQC.
    undef, 1); 
  
  run_system_cmd("mysql --quick -NB -h $options->{host} -u $options->{user} -p$options->{pass} -D $options->{dbname} ".
    '-e "select mf.seq_region_id, mf.seq_region_start, mf.seq_region_end, mf.motif_feature_id, fs.epigenome_id '.
    'FROM motif_feature mf JOIN associated_motif_feature amf USING(motif_feature_id) '.
    'JOIN annotated_feature af USING(annotated_feature_id) JOIN feature_set fs USING(feature_set_id) '.
    "GROUP BY motif_feature_id, epigenome_id ORDER BY seq_region_id, seq_region_start\" > $motifs",
    undef, 1); 

  # Overlap regulatory features with (annotated|motif) features. Store ID pairs into one flat file
  run_system_cmd("bedtools intersect -wa -wb -a $epigenome_regulatory_features -b $annotations ".
    "| awk 'BEGIN {OFS = \"\\t\"} \$5==\$10 {print \$4,\$9,\"annotated\"} ' > $out",
     undef, 1); 

  run_system_cmd("bedtools intersect -wa -wb -a $epigenome_regulatory_features -b $motifs ".
    "| awk 'BEGIN {OFS = \"\\t\"} \$5==\$10 {print \$4,\$9,\"motif\"} ' >> $out",
     undef, 1); 

  # Load into database
  run_system_cmd("mysql -h $options->{host} -u $options->{user} -p$options->{pass} -D $options->{dbname} ".
    '-e "TRUNCATE TABLE regulatory_evidence;"', undef, 1); 

  run_system_cmd("mysql -h $options->{host} -u $options->{user} -p$options->{pass} -D $options->{dbname} ".
    "-e 'LOAD DATA LOCAL INFILE \"$out\" INTO TABLE regulatory_evidence;'", undef, 1);

  unlink $epigenome_regulatory_features;
  unlink $annotations;
  unlink $motifs;
  unlink $out;
}

1;
