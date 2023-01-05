#!/usr/bin/env perl
=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2023] EMBL-European Bioinformatics Institute

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

=head1 create_execution_plan.pl

  create_execution_plan.pl -registry /homes/mnuhn/work_dir_ersa/lib/ensembl-funcgen/registry.pm -species homo_sapiens | less
  create_execution_plan.pl -registry /homes/mnuhn/work_dir_ersa/lib/ensembl-funcgen/registry_new_mouse_encode_data.pm -species mus_musculus | less
  

=head1 SYNOPSIS

=head1 DESCRIPTION

=cut

use strict;
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Logger;
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::ExecutionPlanUtils qw ( 
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

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

$logger->info("registry = " . $registry . "\n");
$logger->info("species  = " . $species  . "\n");

my $ensembl_release_version = '91';
my $root_dir = '';

use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->load_all($registry);

use Bio::EnsEMBL::Funcgen::PeakCallingPlan::ExecutionPlanFactory;

my $execution_plan_factory
  = Bio::EnsEMBL::Funcgen::PeakCallingPlan::ExecutionPlanFactory->new(
    -root_dir                => $root_dir,
    -species                 => $species,
    -ensembl_release_version => $ensembl_release_version,
  );

my $experiment_adaptor = Bio::EnsEMBL::Registry->get_adaptor(
  $species,
  'funcgen',
  'Experiment'
);

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
      and experiment.name not in ("MEL_cell_line_____DNase1_DNase-Seq_ENCODE92")
SQL
;

#       and experiment.name in (
#         "hindbrain_postnatal_0_DNase1_DNase-Seq_ENCODE92",
#         "hindbrain_embryonic_11_5_DNase1_DNase-Seq_ENCODE92",
#         "hindbrain_embryonic_14_5_DNase1_DNase-Seq_ENCODE92",
#         "MEL_cell_line_____DNase1_DNase-Seq_ENCODE92",
#         "midbrain_postnatal_0_DNase1_DNase-Seq_ENCODE92",
#         "liver_embryonic_14_5_DNase1_DNase-Seq_ENCODE92",
#         "neural_tube_embryonic_11_5_DNase1_DNase-Seq_ENCODE92",
#         "stomach_postnatal_0_DNase1_DNase-Seq_ENCODE92",
#         "limb_embryonic_11_5_DNase1_DNase-Seq_ENCODE92",
#         "embryonic_facial_prominence_embryonic_14_5_DNase1_DNase-Seq_ENCODE92",
#         "embryonic_facial_prominence_embryonic_11_5_DNase1_DNase-Seq_ENCODE92",
#         "limb_embryonic_14_5_DNase1_DNase-Seq_ENCODE92",
#         "lung_embryonic_14_5_DNase1_DNase-Seq_ENCODE92",
#         "midbrain_embryonic_14_5_DNase1_DNase-Seq_ENCODE92",
#         "lung_postnatal_0_DNase1_DNase-Seq_ENCODE92",
#         "kidney_postnatal_0_DNase1_DNase-Seq_ENCODE92",
#         "forebrain_postnatal_0_DNase1_DNase-Seq_ENCODE92",
#         "forebrain_embryonic_14_5_DNase1_DNase-Seq_ENCODE92",
#         "heart_postnatal_0_DNase1_DNase-Seq_ENCODE92",
#         "liver_postnatal_0_DNase1_DNase-Seq_ENCODE92"
#       )

# use YAML qw(Dump Bless);
# 
# local $YAML::Indent = 8;
# local $YAML::UseAliases = 0;

$Data::Dumper::Sortkeys = 1;

my @all_execution_plans;

$sql_helper->execute_no_return(
  -SQL          => $do_all_sql,
  -USE_HASHREFS => 1,
  -CALLBACK     => sub {
      my $row = shift;
      
      my $experiment_id = $row->{experiment_id};
      my $experiment = $experiment_adaptor->fetch_by_dbID($experiment_id);
      
      my $execution_plan 
        = $execution_plan_factory->create_execution_plan_for_experiment(
          $experiment
        );

      print Dumper($execution_plan);
      my $execution_plan_expanded = resolve_nonterminal_symbols($execution_plan);
      print summarise($execution_plan_expanded);

      push @all_execution_plans, $execution_plan;
      
      #$execution_plan_adaptor->store($execution_plan_obj);
      return;
    },
);



exit;


my @error_messages;
my %alignment_name_to_experiment_name;
foreach my $execution_plan (@all_execution_plans) {

    my $alignments = $execution_plan->{alignment};
    my @alignment_names = keys %$alignments;
    
    foreach my $alignment_name (@alignment_names) {
    
        my $experiment_name = $alignments->{$alignment_name}->{from_experiment};
    
        if (! exists $alignment_name_to_experiment_name{$alignment_name}) {
            $alignment_name_to_experiment_name{$alignment_name} = $experiment_name;
        }
        if ($alignment_name_to_experiment_name{$alignment_name} ne $experiment_name) {
            push 
                @error_messages, 
                "The alignment name $alignment_name created for the experiment " 
                . $alignment_name_to_experiment_name{$alignment_name} 
                . " is also already being used for the experiment $experiment_name";
        }
    }
}

print "The following errors were detected:\n";
print join "\n", map { "  - " . $_ } @error_messages;

# my @all_alignment_names = keys %alignment_name_to_experiment_name;
# 
# foreach my $alignment_name (@all_alignment_names) {
# 
#     if ($alignment_name_to_experiment_name{$alignment_name} == 1) {
#         delete $alignment_name_to_experiment_name{$alignment_name};
#     }
# }
# 
# print Dumper(\%alignment_name_to_experiment_name);

$logger->finish_log;











