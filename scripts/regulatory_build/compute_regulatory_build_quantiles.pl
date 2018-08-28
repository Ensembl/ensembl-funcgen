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

  A script to compute the number of bases covered by regulatory features both in total and by feature type

  perl scripts/regulatory_build/compute_regulatory_build_quantiles.pl --registry /homes/mnuhn/work_dir_regbuild_script/lib/ensembl-funcgen/registry.pm --species homo_sapiens --tempdir foobar
  perl scripts/regulatory_build/compute_regulatory_build_quantiles.pl --registry /homes/mnuhn/work_dir_regbuild_script/lib/ensembl-funcgen/registry.pm --species mus_musculus --tempdir foobar
  
  perl scripts/regulatory_build/compute_regulatory_build_quantiles.pl --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb6.pm --species homo_sapiens --tempdir foobar

  perl scripts/regulatory_build/compute_regulatory_build_quantiles.pl --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb7.pm --species homo_sapiens --tempdir foobar

=cut

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use Bio::EnsEMBL::Utils::Logger;

my %options;
GetOptions (
    \%options,
    "species|s=s",
    "registry|r=s",
    "tempdir|t=s",
 );

use Hash::Util qw( lock_keys );
lock_keys( %options );

my $species  = $options{'species'};
my $registry = $options{'registry'};

Bio::EnsEMBL::Registry->load_all($registry);
my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');

my $regulatory_build_adaptor = $funcgen_adaptor->get_RegulatoryBuildAdaptor;
my $current_regulatory_build_id = $regulatory_build_adaptor->fetch_current_regulatory_build->dbID;

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

use Statistics::Descriptive;

$logger->info("Computing promoter quantiles\n");
&compute_promoter;

$logger->info("Computing promoter flanking quantiles\n");
&compute_promoter_flanking;

$logger->info("Computing enhancer quantiles\n");
&compute_enhancer;

$logger->info("Computing ctcf quantiles\n");
&compute_ctcf;

$logger->info("Computing transcription factor quantiles\n");
&compute_tf;

$logger->info("Computing open chromatin quantiles\n");
&compute_open_chromatin;

$logger->finish_log;
exit(0);


sub compute_promoter {
  _compute_generic_feature('Promoter', 'promoter');
}
sub compute_promoter_flanking {
  _compute_generic_feature('Promoter Flanking Region', 'promoter_flanking');
}
sub compute_enhancer {
  _compute_generic_feature('Enhancer', 'enhancer');
}
sub compute_ctcf {
  _compute_generic_feature('CTCF Binding Site', 'ctcf');
}
sub compute_tf {
  _compute_generic_feature('TF binding site', 'tf');
}
sub compute_open_chromatin {
  _compute_generic_feature('Open chromatin', 'open_chromatin');
}

sub _compute_generic_feature {

  my $feature_type_name = shift;
  my $statistic_prefix  = shift;

  my $stat = Statistics::Descriptive::Full->new();

  my $sth = $funcgen_adaptor->dbc->db_handle->prepare(
    qq~select (seq_region_end + bound_end_length ) - (seq_region_start-bound_start_length)  + 1 from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "$feature_type_name"~
  );
  $sth->execute;
  my $x = $sth->fetchall_arrayref;

  my $x_flattened = [ map { $_->[0] } @$x ];

  $stat->add_data(@$x_flattened);

  my $skewness = $stat->skewness;
  my $kurtosis = $stat->kurtosis;
  
  if (! defined $skewness) { $skewness = 'null'; }
  if (! defined $kurtosis) { $kurtosis = 'null'; }

  $logger->info("Skewness: " . $skewness . "\n");
  $logger->info("Kurtosis: " . $kurtosis . "\n");

  # https://metacpan.org/pod/Statistics::Descriptive
  # 
  # 0 => zero quartile (Q0) : minimal value
  # 1 => first quartile (Q1) : lower quartile = lowest cut off (25%) of data = 25th percentile
  # 2 => second quartile (Q2) : median = it cuts data set in half = 50th percentile
  # 3 => third quartile (Q3) : upper quartile = highest cut off (25%) of data, or lowest 75% = 75th percentile
  # 4 => fourth quartile (Q4) : maximal value

  my $q0 = $stat->quantile(0);
  my $q1 = $stat->quantile(1);
  my $q2 = $stat->quantile(2);
  my $q3 = $stat->quantile(3);
  my $q4 = $stat->quantile(4);
  
  if (! defined $q0) { $q0 = 0; }
  if (! defined $q1) { $q1 = 0; }
  if (! defined $q2) { $q2 = 0; }
  if (! defined $q3) { $q3 = 0; }
  if (! defined $q4) { $q4 = 0; }

  $logger->info("Q0 " . $q0 . "\n");
  $logger->info("Q1 " . $q1 . "\n");
  $logger->info("Q2 " . $q2 . "\n");
  $logger->info("Q3 " . $q3 . "\n");
  $logger->info("Q4 " . $q4 . "\n");

  my $sql_commands = [
      qq~
      delete from regulatory_build_statistic where statistic in (
        '${statistic_prefix}_q0',
        '${statistic_prefix}_q1',
        '${statistic_prefix}_q2',
        '${statistic_prefix}_q3',
        '${statistic_prefix}_q4',
        '${statistic_prefix}_skewness',
        '${statistic_prefix}_kurtosis'
      );
      ~,
    qq~insert into regulatory_build_statistic (regulatory_build_id, statistic, value) values (${current_regulatory_build_id}, '${statistic_prefix}_q0', $q0);~,
    qq~insert into regulatory_build_statistic (regulatory_build_id, statistic, value) values (${current_regulatory_build_id}, '${statistic_prefix}_q1', $q1);~,
    qq~insert into regulatory_build_statistic (regulatory_build_id, statistic, value) values (${current_regulatory_build_id}, '${statistic_prefix}_q2', $q2);~,
    qq~insert into regulatory_build_statistic (regulatory_build_id, statistic, value) values (${current_regulatory_build_id}, '${statistic_prefix}_q3', $q3);~,
    qq~insert into regulatory_build_statistic (regulatory_build_id, statistic, value) values (${current_regulatory_build_id}, '${statistic_prefix}_q4', $q4);~,

    qq~insert into regulatory_build_statistic (regulatory_build_id, statistic, value) values (${current_regulatory_build_id}, '${statistic_prefix}_skewness', $skewness);~,
    qq~insert into regulatory_build_statistic (regulatory_build_id, statistic, value) values (${current_regulatory_build_id}, '${statistic_prefix}_kurtosis', $kurtosis);~,

  ];
  run_sql_commands($sql_commands);
}

sub run_sql_commands {
  my $sql_commands = shift;

  foreach my $sql_command (@$sql_commands) {
    $logger->info($sql_command);
    $funcgen_adaptor->dbc->do($sql_command);
  }
}
