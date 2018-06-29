#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

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

create_regulatory_evidence.pl

=head1 SYNOPSIS



=head1 DESCRIPTION

=cut

use strict;
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Hash::Util qw( lock_hash );

my $start_time = time;

main();

=head2 main

  Description: Overall process
  Returntype: undef

=cut

sub main {
  print_log("Getting options\n");
  my $options = get_options();

  print_log("Creating regulatory_annotation table\n");
  compute_regulatory_annotations($options);
}

=head2 print_log

  Description: Print conveninence
  Arg1: String
  Returntype: undef
  Side effects: Prints string, with time stamp in front

=cut

sub print_log {
  my ($str) = @_;
  my $runtime = time - $start_time;
  print "[$runtime] $str";
}

=head2 get_options

  Description: Command line options
  Returntype: hashreof

=cut

sub get_options {
  my %options;
  
  use Pod::Usage;
  GetOptions (
    \%options,
    "pass|p=s",
    "port=s",
    "host|h=s",
    "user|u=s",
    "dbname|d=s",
    "dnadb_host=s",
    "dnadb_port=s",
    "dnadb_name=s",
    "dnadb_user=s",
    "dnadb_pass=s",
    "tempdir=s",
  ) or pod2usage( -exitval => 1);
  defined $options{host} || die ("You must define the destination host!\t--host XXXX\n");
  defined $options{user} || die ("You must define the user login!\t--user XXXX\n");
  defined $options{dbname} || die ("You must define the database name!\t--dbname XXXX\n");
  
  lock_hash(%options);
  
  return \%options;
}

sub run {
  my ($cmd) = @_;
  print_log($cmd . "\n");
  use Carp;
  system($cmd) && confess("Failed when running command:\n$cmd\n");
}

=head2 compute_regulatory_annotations

  Description: Assign motifs and annotations to regulatory features
  Arg1: options hash ref
  Returntype: undef
  Side effects: writes into regulatory_evidences table

=cut

sub compute_regulatory_annotations {

  my $options = shift;

  my $tempdir = $options->{tempdir};
  
  use File::Path qw( mkpath );
  mkpath ($tempdir);
  
  my $regulatory_features        = $tempdir . '/' . 'regulatory_features.bed';
  my $peaks                      = $tempdir . '/' . 'peaks.bed';
  my $motifs                     = $tempdir . '/' . 'motifs.bed';
  my $regulatory_evidence_peaks  = $tempdir . '/' . 'regulatory_evidence_peaks.bed';
  my $regulatory_evidence_motifs = $tempdir . '/' . 'regulatory_evidence_motifs.bed';
  
  # Extract regulatory, annotated, and motif features into flat tab-delimited files.
  # The common format of these files is:
  # chrom	start	end	ID
  
  # "convert" in the "order by" clause makes MySql sort lexicographically. 
  # This is necessary, or bedtools will complain later.
  #
  
  run(
    qq(
	mysql --quick -NB -h $options->{host} -u $options->{user} -p$options->{pass} -P $options->{port} -D $options->{dbname} -e '
	  select
	    rf.seq_region_id,
	    rf.seq_region_start - rf.bound_start_length as start,
	    rf.seq_region_end + rf.bound_end_length as end,
	    regulatory_activity_id,
	    epigenome_id
	  from
	    regulatory_feature rf
	    join regulatory_activity using (regulatory_feature_id)
	  order by
	    convert(rf.seq_region_id, char),
	    start;
	' \\
 	> $regulatory_features
    )
  );

  run(
    qq(
	mysql --quick -NB -h $options->{host} -u $options->{user} -p$options->{pass} -P $options->{port} -D $options->{dbname} -e '
      select
        peak.seq_region_id,
        peak.seq_region_start,
        peak.seq_region_end,
        peak.peak_id,
        peak_calling.epigenome_id
      FROM
        peak
        join peak_calling USING (peak_calling_id)
      ORDER BY
        convert(peak.seq_region_id,    char),
        peak.seq_region_start;
	' \\
	> $peaks
    )
  );

  # run(
  #   qq(
	# mysql --quick -NB -h $options->{host} -u $options->{user} -p$options->{pass} -D $options->{dbname} -e '
	#   select
	#     mf.seq_region_id, mf.seq_region_start, mf.seq_region_end, mf.motif_feature_id, fs.epigenome_id
	#   from
	#     motif_feature mf
	#     JOIN associated_motif_feature amf USING(motif_feature_id)
	#     JOIN annotated_feature af USING(annotated_feature_id)
	#     JOIN feature_set fs USING(feature_set_id)
	#   GROUP BY
	#     motif_feature_id,
	#     epigenome_id
	#   ORDER BY
	#     seq_region_id,
	#     seq_region_start;
	#  ' > $motifs
  #   )
  # );

  run(
    qq(
    mysql --quick -NB -h $options->{host} -u $options->{user} -p$options->{pass} -P $options->{port} -D $options->{dbname} -e '
      select
        mf.seq_region_id, mf.seq_region_start, mf.seq_region_end, mf.motif_feature_id, "PHONY" 
      from
        motif_feature mf
      ORDER BY
        convert(seq_region_id,    char),
        seq_region_start;
  ' > $motifs
    )
  );

  # Overlap regulatory features with (annotated|motif) features. Store ID pairs into one flat file
  run("bedtools intersect -sorted -wa -wb -a $regulatory_features -b $peaks  | awk 'BEGIN {OFS = \"\\t\"} \$5==\$10 {print \$4,\$9,\"annotated\"} ' > $regulatory_evidence_peaks");
  run("bedtools intersect -sorted -wa -wb -a $regulatory_features -b $motifs | awk 'BEGIN {OFS = \"\\t\"} {print \$4,\$9,\"motif\"} ' > $regulatory_evidence_motifs");

  # Load into database
  run("mysql -h $options->{host} -u $options->{user} -p$options->{pass} -D $options->{dbname} -P $options->{port} -e 'TRUNCATE TABLE regulatory_evidence;'");
  run("mysql -h $options->{host} -u $options->{user} -p$options->{pass} -D $options->{dbname} -P $options->{port} -e 'LOAD DATA LOCAL INFILE \"$regulatory_evidence_peaks\"  INTO TABLE regulatory_evidence;'");
  run("mysql -h $options->{host} -u $options->{user} -p$options->{pass} -D $options->{dbname} -P $options->{port} -e 'LOAD DATA LOCAL INFILE \"$regulatory_evidence_motifs\" INTO TABLE regulatory_evidence;'");

}

