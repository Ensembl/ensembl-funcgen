#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

$logger->info("Computing total regulatory features counts\n");
compute_total_counts($funcgen_adaptor);

$logger->info("Computing averages\n");
compute_averages($funcgen_adaptor);

$logger->info("Computing lengths\n");
compute_lengths($funcgen_adaptor);

$logger->finish_log;
exit(0);

sub run_sql_commands {
  my $sql_commands = shift;

  foreach my $sql_command (@$sql_commands) {
    $logger->info($sql_command);
    $funcgen_adaptor->dbc->do($sql_command);
  }
}

sub compute_lengths {

  my $funcgen_adaptor = shift;

  my $sql_commands = [
    qq~
    delete from regulatory_build_statistic where statistic in (
      'sum_length_promoter',
      'sum_length_enhancer',
      'sum_length_promoter_flanking_region',
      'sum_length_transcription_factor_binding_site',
      'sum_length_open_chromatin',
      'sum_length_ctcf_binding_site'
    );
    ~,
    qq~
    insert into regulatory_build_statistic (regulatory_build_id, statistic, value) (
        select regulatory_build_id, 'sum_length_promoter', SUM( (seq_region_end + bound_end_length ) - (seq_region_start-bound_start_length)  + 1) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "Promoter" group by feature_type.name, regulatory_build_id
    );
    ~,
    qq~
    insert into regulatory_build_statistic (regulatory_build_id, statistic, value) (
        select regulatory_build_id, 'sum_length_enhancer', SUM( (seq_region_end + bound_end_length ) - (seq_region_start-bound_start_length)  + 1) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "Enhancer" group by feature_type.name, regulatory_build_id
    );
    ~,
    qq~
    insert into regulatory_build_statistic (regulatory_build_id, statistic, value) (
        select regulatory_build_id, 'sum_length_promoter_flanking_region', SUM( (seq_region_end + bound_end_length ) - (seq_region_start-bound_start_length)  + 1) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "Promoter Flanking Region" group by feature_type.name, regulatory_build_id
    );
    ~,
    qq~
    insert into regulatory_build_statistic (regulatory_build_id, statistic, value) (
        select regulatory_build_id, 'sum_length_transcription_factor_binding_site', SUM( (seq_region_end + bound_end_length ) - (seq_region_start-bound_start_length)  + 1) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "TF binding site" group by feature_type.name, regulatory_build_id
    );
    ~,
    qq~
    insert into regulatory_build_statistic (regulatory_build_id, statistic, value) (
        select regulatory_build_id, 'sum_length_open_chromatin', SUM( (seq_region_end + bound_end_length ) - (seq_region_start-bound_start_length)  + 1) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "Open chromatin" group by feature_type.name, regulatory_build_id
    );
    ~,
    qq~
    insert into regulatory_build_statistic (regulatory_build_id, statistic, value) (
        select regulatory_build_id, 'sum_length_ctcf_binding_site', SUM( (seq_region_end + bound_end_length ) - (seq_region_start-bound_start_length)  + 1) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "CTCF Binding Site" group by feature_type.name, regulatory_build_id
    );
    ~

    ];
  run_sql_commands($sql_commands)
}

sub compute_total_counts {

  my $funcgen_adaptor = shift;

  my $sql_commands = [
    qq~
    delete from regulatory_build_statistic where statistic in (
      'number_regulatory_features',
      'number_promoter',
      'number_enhancer',
      'number_promoter_flanking_region',
      'number_transcription_factor_binding_site',
      'number_open_chromatin',
      'number_ctcf_binding_site'
    );
    ~,
    "
    insert into regulatory_build_statistic (regulatory_build_id, statistic, value) (
        select 
        regulatory_build_id, 
        'number_regulatory_features', 
        count(regulatory_feature_id)
        from 
        regulatory_feature 
        group by regulatory_build_id
    );
    ",
    qq~
    insert into regulatory_build_statistic (regulatory_build_id, statistic, value) (
        select regulatory_build_id, 'number_promoter', count(regulatory_feature_id) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "Promoter" group by feature_type.name, regulatory_build_id
    );
    ~,
    qq~
    insert into regulatory_build_statistic (regulatory_build_id, statistic, value) (
        select regulatory_build_id, 'number_enhancer', count(regulatory_feature_id) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "Enhancer" group by feature_type.name, regulatory_build_id
    );
    ~,
    qq~
    insert into regulatory_build_statistic (regulatory_build_id, statistic, value) (
        select regulatory_build_id, 'number_promoter_flanking_region', count(regulatory_feature_id) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "Promoter Flanking Region" group by feature_type.name, regulatory_build_id
    );
    ~,
    qq~
    insert into regulatory_build_statistic (regulatory_build_id, statistic, value) (
        select regulatory_build_id, 'number_transcription_factor_binding_site', count(regulatory_feature_id) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "TF binding site" group by feature_type.name, regulatory_build_id
    );
    ~,
    qq~
    insert into regulatory_build_statistic (regulatory_build_id, statistic, value) (
        select regulatory_build_id, 'number_open_chromatin', count(regulatory_feature_id) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "Open chromatin" group by feature_type.name, regulatory_build_id
    );
    ~,
    qq~
    insert into regulatory_build_statistic (regulatory_build_id, statistic, value) (
        select regulatory_build_id, 'number_ctcf_binding_site', count(regulatory_feature_id) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "CTCF Binding Site" group by feature_type.name, regulatory_build_id
    );
    ~,
    ];
  run_sql_commands($sql_commands)
}

sub compute_averages {

  my $funcgen_adaptor = shift;

  my $sql_commands = [
    qq~
    delete from regulatory_build_statistic where statistic in (
      'average_length_promoter',
      'average_length_enhancer',
      'average_length_promoter_flanking_region',
      'average_length_transcription_factor_binding_site',
      'average_length_open_chromatin',
      'average_length_ctcf_binding_site'
    );
    ~,
    qq~
    insert into regulatory_build_statistic (regulatory_build_id, statistic, value) (
        select regulatory_build_id, 'average_length_promoter', AVG( (seq_region_end + bound_end_length ) - (seq_region_start-bound_start_length)  + 1) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "Promoter" group by feature_type.name, regulatory_build_id
    );
    ~,
    qq~
    insert into regulatory_build_statistic (regulatory_build_id, statistic, value) (
        select regulatory_build_id, 'average_length_enhancer', AVG( (seq_region_end + bound_end_length ) - (seq_region_start-bound_start_length)  + 1) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "Enhancer" group by feature_type.name, regulatory_build_id
    );
    ~,
    qq~
    insert into regulatory_build_statistic (regulatory_build_id, statistic, value) (
        select regulatory_build_id, 'average_length_promoter_flanking_region', AVG( (seq_region_end + bound_end_length ) - (seq_region_start-bound_start_length)  + 1) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "Promoter Flanking Region" group by feature_type.name, regulatory_build_id
    );
    ~,
    qq~
    insert into regulatory_build_statistic (regulatory_build_id, statistic, value) (
        select regulatory_build_id, 'average_length_transcription_factor_binding_site', AVG( (seq_region_end + bound_end_length ) - (seq_region_start-bound_start_length)  + 1) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "TF binding site" group by feature_type.name, regulatory_build_id
    );
    ~,
    qq~
    insert into regulatory_build_statistic (regulatory_build_id, statistic, value) (
        select regulatory_build_id, 'average_length_open_chromatin', AVG( (seq_region_end + bound_end_length ) - (seq_region_start-bound_start_length)  + 1) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "Open chromatin" group by feature_type.name, regulatory_build_id
    );
    ~,
    qq~
    insert into regulatory_build_statistic (regulatory_build_id, statistic, value) (
        select regulatory_build_id, 'average_length_ctcf_binding_site', AVG( (seq_region_end + bound_end_length ) - (seq_region_start-bound_start_length)  + 1) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "CTCF Binding Site" group by feature_type.name, regulatory_build_id
    );
    ~,
  ];
  run_sql_commands($sql_commands)
}

