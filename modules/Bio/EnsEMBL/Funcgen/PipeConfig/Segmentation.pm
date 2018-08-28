=pod 

=head1 NAME

    Bio::EnsEMBL::Funcgen::PipeConfig::Segmentation

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 LICENSE

    Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
    Copyright [2016-2017] EMBL-European Bioinformatics Institute

    Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

         http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software distributed under the License
    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and limitations under the License.

=head1 CONTACT

=cut
package Bio::EnsEMBL::Funcgen::PipeConfig::Segmentation;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base 'Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::ResourceClasses';
#use base ('Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::Base');

sub pipeline_wide_parameters {
  my $self = shift;
  
  return {
    %{$self->SUPER::pipeline_wide_parameters},
    pipeline_name            => $self->o('pipeline_name'),
    tempdir                  => $self->o('tempdir'),
    tempdir_segmentation     => $self->o('tempdir') . '/segmentation',
    tempdir_regulatory_build => $self->o('tempdir') . '/regulatory_build',
    data_root_dir            => $self->o('data_root_dir'),
    reference_data_root_dir  => $self->o('reference_data_root_dir'),
    reg_conf                 => $self->o('reg_conf'),
    #ChromHMM                 => '/nfs/production/panda/ensembl/funcgen/ChromHMM/ChromHMM.jar',
    ChromHMM                 => '/nfs/production/panda/ensembl/funcgen/ChromHMM/1.17/ChromHMM/ChromHMM.jar',
  };
}

sub pipeline_analyses {
    my ($self) = @_;
    return [
      {   -logic_name => 'start_segmentation',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -flow_into => { 
            MAIN => 'truncate_regulatory_build_tables',
          },
      },
      {
          -logic_name  => 'truncate_regulatory_build_tables',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
          -parameters => {
              sql     => [
                'truncate regulatory_build;',
                'truncate regulatory_feature;',
                'truncate regulatory_activity;',
                'truncate regulatory_build_epigenome;',
                'truncate regulatory_evidence;',
                'truncate regulatory_build_statistic;',
                'truncate segmentation_file;',
                'delete from data_file where table_name = "segmentation_file";',
                'truncate segmentation_state_assignment;',
                'truncate segmentation_state_emission;',
                'truncate segmentation_cell_tables',
                'truncate segmentation',
              ],
              db_conn => 'funcgen:#species#',
          },
          -flow_into   => {
              MAIN => 'segmentation_job_factory',
          },
      },
      {   -logic_name => 'segmentation_job_factory',
          -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::Segmentation::SegmentationJobFactory',
          -parameters => {
            tempdir => '#tempdir_segmentation#/#species#/',
          },
          -flow_into => {
            2 => 'create_cell_tables',
          },
      },
      {   -logic_name => 'create_cell_tables',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
              cmd => 
                  q( 
                    classify_epigenome_to_segmentation_run.pl \
                        --registry #reg_conf# \
                        --species #species# \
                        --partition_by_experimental_group 1
                  )
          },
          -flow_into => { 
            MAIN => 'seed_segmentation_jobs',
          },
      },
      {
          -logic_name  => 'seed_segmentation_jobs',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
          -parameters => {
              db_conn    => 'funcgen:#species#',
              inputquery => '
                select 
                    distinct superclass, 
                    class
                from 
                    segmentation_cell_tables
              ',
          },
          -flow_into   => {
            '2->A'     => { 
                'make_segmentation_dir' => INPUT_PLUS({
                    class                 => '#class#',
                    superclass            => '#superclass#',
                    segmentation_name     => '#superclass#_#class#',
                    celltable_file        => 'celltable.#superclass#.#class#.txt',
                    learn_model_directory => '#learnmodel_output_dir#/#superclass#/#class#/',
                })
            },
            'A->MAIN'  => { 
                'generate_segmentation_parameter_file'  => INPUT_PLUS({
                    foo => 'bar',
                })
            },
          },
      },
      {   -logic_name => 'make_segmentation_dir',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          #-module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -parameters => {
              cmd => 
                  q( rm    -rf #tempdir_segmentation#/#species#/#superclass#/#class# ; )
                . q( mkdir -p  #tempdir_segmentation#/#species#/#superclass#/#class#   ),
          },
          -flow_into => { 
            MAIN => 'generate_cell_table_file',
          },
      },
      {   -logic_name => 'generate_cell_table_file',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
              cmd => 
                  q( 
                    generate_cell_table_file.pl        \
                        --registry        #reg_conf#   \
                        --species         #species#    \
                        --superclass      #superclass# \
                        --class           #class#      \
                        --cell_table_file #tempdir_segmentation#/#species#/#superclass#/#class#/celltable.#superclass#.#class#.txt
                  )
          },
          -flow_into => { 
            MAIN => [ 
                'binarize', 
                #'record_segmentation_as_done' 
            ],
          },
      },
      {   -logic_name => 'binarize',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -max_retry_count => 0,
          -rc_name    => 'binarization',
          -parameters => { 
            cmd => q( java -Xmx14500m                     )
                 . q(   -jar #ChromHMM#                   )
                 . q(   BinarizeBam                       )
                 . q(   #chromosome_length_file#          )
                 . q(   #data_root_dir_species_assembly#  )
                 . q(   #tempdir_segmentation#/#species#/#superclass#/#class#/celltable.#superclass#.#class#.txt )
                 . q(   #binarized_bam_dir#/#superclass#/#class#/ ),
          },
          -flow_into => {
            MAIN     => 'delete_non_chromosomal_files_from_binarization_dir',
          },
      },
      {   -logic_name => 'delete_non_chromosomal_files_from_binarization_dir',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => { 
            cmd => 
                   q( 
                    delete_non_chromosomal_files_from_binarization_dir.pl \
                        --registry  #reg_conf#   \
                        --species   #species#    \
                        --directory #binarized_bam_dir#/#superclass#/#class#
                  ),
          },
          -flow_into => {
            MAIN     => 'learn_model',
          },
      },
      {   -logic_name => 'learn_model',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -max_retry_count => 0,
          -rc_name    => 'learn_model',
          -parameters => {
            cmd     => 'run_ChromHMM_skip_nonsense_errors.pl --cmd "#run_cmd#"',
            run_cmd => 
                q( export DISPLAY= ;             )
                . q( java -Xmx30000m             )
                . q( -jar #ChromHMM#             )
                . q( LearnModel                  )
                . q( -r 300 -p 12                )
                . q( -l #chromosome_length_file# )
                . q( #binarized_bam_dir#/#superclass#/#class#/         )
                . q( #learn_model_directory#     )
                . q( 25                          )
                . q( #assembly#                  )
          },
          -flow_into => {
            MAIN     => 'record_segmentation_as_done',
          },
      },
      {   -logic_name => 'record_segmentation_as_done',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -parameters => {
            record => {
                learn_model_directory => '#learn_model_directory#',
                segmentation_name     => '#superclass#_#class#',
                celltable_file        => 'celltable.#superclass#.#class#.txt',
            }
          },
          -flow_into => {
            MAIN => '?accu_name=segmentation_lists&accu_address={superclass}{class}&accu_input_variable=record'
          },
      },
      {   -logic_name => 'generate_segmentation_parameter_file',
          -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::Segmentation::GenerateSegmentationParameterFile',
          -parameters => {
            file => '#segmentation_parameter_file#',
          },
          -flow_into => {
            MAIN     => { 
                'make_regbuild_dir' => INPUT_PLUS({ 
                    segmentation_lists => '#segmentation_lists#' 
                 })
            },
          },
      },
      {   -logic_name => 'make_regbuild_dir',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
              cmd => 
                  q(rm    -rf #tempdir_regulatory_build#/#species#/#assembly# ; )
                . q(mkdir -p  #tempdir_regulatory_build#/#species#/#assembly#   ),
          },
          -flow_into => { 
            MAIN     => 'build_regulatory_features',
          },
      },
      {   -logic_name => 'build_regulatory_features',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -rc_name    => '4Gb_job',
          -parameters => {
            cmd => 
                q( new_build_regulatory_features.pl             )
                . q( --species #species#                        )
                . q( --out #tempdir_regulatory_build#/#species# )
                . q( -dump #segmentation_parameter_file#     )
                . q( -a #assembly#                )
                . q( --db   #funcgen_dbname#      )
                . q( --user #funcgen_username#    )
                . q( --host #funcgen_host#        )
                . q( --pass #funcgen_password#    )
                . q( --port #funcgen_port#        )
                . q( --dnadb_name #core_dbname#   )
                . q( --dnadb_host #core_host#     )
                . q( --dnadb_user #core_username# )
                . q( --dnadb_pass #core_password# )
                . q( --dnadb_port #core_port#     )
                . q( --chrom_lengths #chromosome_length_file# )
          },
          -flow_into => {
            MAIN => [ 
                'load_build', 
                'seed_load_assignment_jobs',
            ]
          },
      },
      {   -logic_name => 'load_build',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -rc_name    => '32Gb_job_2cpus',
          -parameters => {
            cmd => 
                q( load_build.pl )
                . q( --dbname   #funcgen_dbname# )
                . q( --user     #funcgen_username#    )
                . q( --host     #funcgen_host#        )
                . q( --pass     #funcgen_password#    )
                . q( --port     #funcgen_port#        )
                . q( --dnadb_name #core_dbname#   )
                . q( --dnadb_host #core_host#     )
                . q( --dnadb_user #core_username# )
                . q( --dnadb_pass #core_password# )
                . q( --dnadb_port #core_port#     )
                . q( --base_dir #tempdir_regulatory_build#/#species#/#assembly# )
          },
          -flow_into => {
            MAIN => 'populate_meta_coord'
          },
      },
      {   -logic_name => 'populate_meta_coord',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          #-module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -parameters => {
            cmd => qq( populate_meta_coord.pl    )
              . qq( --species  #species#         )
              . qq( --registry #reg_conf#        )
          },
          -flow_into => {
            MAIN => 'populate_regulatory_build_epigenome_table',
          },
      },
      {
          -logic_name  => 'seed_load_assignment_jobs',
          -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::Segmentation::SeedLoadAssignmentJobs',
          -parameters => {
              db_conn    => 'funcgen:#species#',
          },
          -flow_into   => {
            '2->A'     => 'load_state_emissions',
            'A->MAIN'  => 'load_state_assignments',
          },
      },
      {   -logic_name => 'load_state_emissions',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 
                q(
                  load_state_emissions.pl \
                      --species          #species# \
                      --registry         #reg_conf# \
                      --emissions_file   #learn_model_directory#/emissions_25.txt \
                      --segmentation     #segmentation_name#
                )
          },
      },
      {   -logic_name => 'load_state_assignments',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 
                q(
                  load_state_assignments.pl \
                      --species          #species# \
                      --registry         #reg_conf# \
                      --assignments_file #tempdir_regulatory_build#/#species#/tmp/assignments.txt
                )
          },
      },
      {
          -logic_name  => 'populate_regulatory_build_epigenome_table',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
          -parameters => {
              sql     => [
                "truncate regulatory_build_epigenome;",
                "insert into regulatory_build_epigenome (
                    regulatory_build_id, 
                    epigenome_id
                 ) select 
                      distinct regulatory_build_id, epigenome_id 
                   from 
                    regulatory_build 
                    join regulatory_feature using (regulatory_build_id) 
                    join regulatory_activity using (regulatory_feature_id);
                ",
              ],
              db_conn => 'funcgen:#species#',
          },
          -flow_into => {
            MAIN => 'seed_register_segmentation_files_jobs',
          },
      },
      {
          -logic_name  => 'seed_register_segmentation_files_jobs',
          -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::Segmentation::SeedLoadAssignmentJobs',
          -parameters => {
              db_conn    => 'funcgen:#species#',
          },
          -flow_into   => {
            2 => 'register_segmentation_files',
          },
      },
      {   -logic_name => 'register_segmentation_files',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 
                q(
                  register_segmentation_files.pl  \
                      --species          #species# \
                      --registry         #reg_conf# \
                      --segmentation_directory #tempdir_regulatory_build#/#species#/#assembly#/segmentations/#segmentation_name# \
                      --db_file_species_assembly_dir     #data_root_dir#/#species#/#assembly#/funcgen/segmentation_file/#ensembl_release_version# \
                      --db_file_relative_dir             /funcgen/segmentation_file/#ensembl_release_version# \
                      --segmentation_name                #segmentation_name#
                )
          },
      },
    ];
}

1;
