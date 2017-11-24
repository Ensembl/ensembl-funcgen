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
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either 405ress or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 create_execution_plan_new.pl

  create_execution_plan_new.pl -registry /homes/mnuhn/work_dir_ersa/lib/ensembl-funcgen/registry.pm -species homo_sapiens

=head1 SYNOPSIS

=head1 DESCRIPTION

=cut

use strict;
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Logger;

use Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::SignalFilePlanBuilder;
use Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::RemoveDuplicatesPlanFactory;
use Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::ExecutionPlanUtils qw ( 
  resolve_nonterminal_symbols
  summarise
);

my $registry;
my $species;

my %config_hash = (
  "registry" => \$registry,
  "species"  => \$species,
);

my $result = GetOptions(
  \%config_hash,
  'registry=s',
  'species=s',
);

$Data::Dumper::Deepcopy = 1;

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

$logger->info("registry = " . $registry . "\n");
$logger->info("species  = " . $species  . "\n");

my $ensembl_release_version = '91';

$logger->info("ensembl_release_version = " . $ensembl_release_version . "\n");

my $root_dir = '';

use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->load_all($registry);

my $funcgen_db_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
my $core_db_adaptor    = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'core');

my $experiment_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'Experiment');
my $coordsystem_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'coordsystem');

# my $regulatory_build_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'RegulatoryBuild');
# my $regulatory_build = $regulatory_build_adaptor->fetch_current_regulatory_build;
# 
# print Dumper($regulatory_build);
# 
# exit;

my $default_chromosome_coordsystem = $coordsystem_adaptor->fetch_by_name('chromosome');
my $default_assembly = $default_chromosome_coordsystem->version;

$logger->info("assembly = " . $default_assembly . "\n");

use Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::DirectoryNameBuilder;
my $directory_name_builder 
  = Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::DirectoryNameBuilder
    ->new(
      -root_dir                => $root_dir,
      -species                 => $species,
      -assembly                => $default_assembly,
      -ensembl_release_version => $ensembl_release_version,
    );

use Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::Director;
my $chip_seq_analysis_director = Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::Director->new;

use Bio::EnsEMBL::Utils::SqlHelper;
my $sql_helper = Bio::EnsEMBL::Utils::SqlHelper->new(
  -DB_CONNECTION => $experiment_adaptor->dbc
);

my $do_all_sql =<<SQL
    select 
      experiment_id 
    from 
      experiment 
      join feature_type using (feature_type_id) 
      join experimental_group using (experimental_group_id) 
    where 
      is_control = 0 
      and class in (
        "Histone", 
        "Open Chromatin",
        "Transcription Factor",
        "Polymerase"
      )
      and experimental_group.name != "BLUEPRINT"
SQL
;

my $do_one_idr_type_each_sql =<<SQL
    select 
      experiment_id 
    from 
      experiment 
      join feature_type using (feature_type_id) 
    where 
      experiment.name in (
        'iPS_20b_H3K27me3_ChIP-Seq_Roadmap85',
        'H1_neuronal_progenitor_H3K27ac_ChIP-Seq_Roadmap85',
        'Psoas_Muscle_H3K4me3_ChIP-Seq_Roadmap85'
      )
SQL
;

use Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::Director;
my $chip_seq_analysis_director = Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::Director->new;

$sql_helper->execute_no_return(
  -SQL          => $do_all_sql,
  -USE_HASHREFS => 1,
  -CALLBACK     => sub {
      my $row = shift;
      
      my $experiment_id = $row->{experiment_id};
      my $experiment = $experiment_adaptor->fetch_by_dbID($experiment_id);
      
      my $execution_plan 
        = $chip_seq_analysis_director->construct_execution_plan(
          {
            species                => $species, 
            assembly               => $default_assembly, 
            experiment             => $experiment,
            directory_name_builder => $directory_name_builder
          }
        );
        
      print summarise($execution_plan);
      my $execution_plan_expanded = resolve_nonterminal_symbols($execution_plan);
      print summarise($execution_plan_expanded);
      return;
    },
);

$logger->finish_log;

