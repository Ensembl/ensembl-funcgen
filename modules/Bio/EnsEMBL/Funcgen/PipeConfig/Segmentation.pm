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
use base ('Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::Base');

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
    ChromHMM                 => '/nfs/production/panda/ensembl/funcgen/ChromHMM/ChromHMM.jar',
  };
}

sub pipeline_analyses {
    my ($self) = @_;
    return [
      {   -logic_name => 'start_segmentation',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -flow_into => { 
            'MAIN->A' => 'truncate_regulatory_build_tables',
            'A->MAIN' => 'segmentation_done'
          },
      },
      {
          -logic_name  => 'truncate_regulatory_build_tables',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
          -parameters => {
              sql     => [
                "truncate regulatory_build;",
                "truncate regulatory_feature;",
                "truncate regulatory_activity;",
                "truncate regulatory_build_epigenome;",
                "truncate regulatory_evidence;",
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
            2 => 'populate_meta_coord',
          },
      },
      {   -logic_name => 'populate_meta_coord',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => qq( populate_meta_coord.pl    )
              . qq( --species  #species#         )
              . qq( --registry #reg_conf#        )
          },
          -flow_into => {
            MAIN     => 'generate_cell_table_file',
          },
      },
      {   -logic_name => 'generate_cell_table_file',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => qq( generate_cell_table_file.pl      )
              . qq( --species         #species#         )
              . qq( --registry        #reg_conf#        )
              . qq( --cell_table_file #tempdir_segmentation#/#species#/celltable.txt )
          },
          -flow_into => {
            MAIN     => 'binarize',
          },
      },
      {   -logic_name => 'binarize',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -rc_name    => 'binarization',
          -parameters => { 
            cmd => qq(java -Xmx14500m )
              . qq( -jar #ChromHMM# )
              . qq( BinarizeBam                       )
              . qq( #chromosome_length_file#          )
              . qq( #data_root_dir_species_assembly#  )
              . qq( #tempdir_segmentation#/#species#/celltable.txt )
              . qq( #binarized_bam_dir#               ),
          },
          -flow_into => {
            MAIN     => 'generate_segmentation_parameter_file',
          },
      },
      {   -logic_name => 'generate_segmentation_parameter_file',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => q( echo -e "segmentation\t#species#_segmentation\tChromHMM\t#learnmodel_output_dir#" > #segmentation_parameter_file# )
          },
          -flow_into => { 
            MAIN     => 'learn_model',
          },
      },
      {   -logic_name => 'learn_model',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -rc_name    => 'learn_model',
          -parameters => {
            cmd => 
                qq( export DISPLAY= ;             )
                . qq( java -Xmx30000m             )
                . qq( -jar #ChromHMM#             )
                . qq( LearnModel                  )
                . qq( -r 300 -p 12                )
                . qq( -l #chromosome_length_file# )
                . qq( #binarized_bam_dir#         )
                . qq( #learnmodel_output_dir#     )
                . qq( 25                          )
                . qq( #assembly#                  )
          },
          -flow_into => {
            MAIN => 'make_regbuild_dir',
          },
      },
      {   -logic_name => 'make_regbuild_dir',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
              cmd => qq(mkdir -p #tempdir_regulatory_build#/#species#/#assembly#),
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
                qq( new_build_regulatory_features.pl             )
                . qq( --species #species#                        )
                . qq( --out #tempdir_regulatory_build#/#species# )
                . qq( -dump #segmentation_parameter_file#     )
                . qq( -a #assembly#                )
                . qq( --db   #funcgen_dbname#      )
                . qq( --user #funcgen_username#    )
                . qq( --host #funcgen_host#        )
                . qq( --pass #funcgen_password#    )
                . qq( --port #funcgen_port#        )
                . qq( --dnadb_name #core_dbname#   )
                . qq( --dnadb_host #core_host#     )
                . qq( --dnadb_user #core_username# )
                . qq( --dnadb_pass #core_password# )
                . qq( --dnadb_port #core_port#     )
                . qq( --chrom_lengths #chromosome_length_file# )
          },
          -flow_into => {
            MAIN => 'load_build',
          },
      },
      {   -logic_name => 'load_build',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -rc_name    => '4Gb_job_2cpus',
          -parameters => {
            cmd => 
                qq( load_build.pl )
                . qq( --dbname   #funcgen_dbname# )
                . qq( --user     #funcgen_username#    )
                . qq( --host     #funcgen_host#        )
                . qq( --pass     #funcgen_password#    )
                . qq( --port     #funcgen_port#        )
                . qq( --dnadb_name #core_dbname#   )
                . qq( --dnadb_host #core_host#     )
                . qq( --dnadb_user #core_username# )
                . qq( --dnadb_pass #core_password# )
                . qq( --dnadb_port #core_port#     )
                . qq( --base_dir #tempdir_regulatory_build#/#species#/#assembly# )
          },
          -flow_into => {
            MAIN => 'create_regulatory_build_statistics',
          },
      },
      {
          -logic_name  => 'create_regulatory_build_statistics',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
          -parameters => {
              sql     => [
                "drop table if exists regulatory_build_statistics",

                "CREATE TABLE regulatory_build_statistics (
                  regulatory_build_statistics_id INT(30) UNSIGNED NOT NULL AUTO_INCREMENT,
                  regulatory_build_id INT(22) UNSIGNED NOT NULL,
                  statistic                VARCHAR(128) NOT NULL,
                  value                    BIGINT(11) UNSIGNED DEFAULT '0' NOT NULL,
                  PRIMARY KEY (regulatory_build_statistics_id),
                  UNIQUE KEY stats_uniq(statistic, regulatory_build_id)
                )",

                "insert into regulatory_build_statistics (regulatory_build_id, statistic, value) (
                  select 
                    regulatory_build_id, 
                    'Number of regulatory features', 
                    count(regulatory_feature_id)
                  from 
                    regulatory_feature 
                  group by regulatory_build_id
                )",

                "insert into regulatory_build_statistics (regulatory_build_id, statistic, value) (
                  select 
                    regulatory_build_id,
                    feature_type.name, 
                    count(regulatory_feature_id)
                  from 
                    regulatory_feature 
                    join feature_type using (feature_type_id) 
                  group by 
                    feature_type.name, 
                    regulatory_build_id
                )",
              ],
              db_conn => 'funcgen:#species#',
          },
          -flow_into => {
            MAIN => 'segmentation_done',
          },
      },
#       {   -logic_name => 'create_regulatory_evidence',
#           -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
#           -parameters => {
#             cmd => 
#                 qq( create_regulatory_evidence.pl )
#                 . qq( --dbname   #funcgen_dbname# )
#                 . qq( --user     #funcgen_username#    )
#                 . qq( --host     #funcgen_host#        )
#                 . qq( --pass     #funcgen_password#    )
#                 . qq( --port     #funcgen_port#        )
#                 . qq( --dnadb_name #core_dbname#   )
#                 . qq( --dnadb_host #core_host#     )
#                 . qq( --dnadb_user #core_username# )
#                 . qq( --dnadb_pass #core_password# )
#                 . qq( --dnadb_port #core_port#     )
#                 . qq( --tempdir #tempdir_regulatory_build#/#species#/#assembly# )
#           },
#       },
      {   -logic_name => 'segmentation_done',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      },
    ];
}

sub resource_classes {
    my ($self) = @_;
    return {
         %{$self->SUPER::resource_classes},  # inherit 'default' from the parent class
         'binarization' => {'LSF' => '-M16000 -R"select[mem>16000] rusage[mem=16000]" -n 5'                },
         'learn_model'  => {'LSF' => '-M31000 -R"span[hosts=1] select[mem>31000] rusage[mem=31000]" -n 12' },
    };
}

1;
