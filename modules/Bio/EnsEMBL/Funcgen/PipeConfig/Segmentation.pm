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
            MAIN => 'truncate_regulatory_build_tables',
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
                "truncate regulatory_build_statistic;",
                "truncate segmentation_file;",
                "delete from data_file where table_name = 'segmentation_file';",
                "truncate segmentation_state_assignment;",
                "truncate segmentation_state_emission;",
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
            MAIN     => 'create_cell_tables',
          },
      },
      {
          -logic_name  => 'create_cell_tables',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
          -parameters => {
              sql     => [
                q~ drop table if exists temp_segmentation_candidates ~,
                q~
                    create table temp_segmentation_candidates as 
                    select
                        epigenome.production_name as epigenome,
                        feature_type.name as feature_type,
                        signal_bam.path as signal_bam_path,
                        control_bam.path as control_bam_path
                    from
                        experiment
                        join epigenome using (epigenome_id)
                        join feature_type using (feature_type_id)
                        join alignment signal_alignment on (
                        experiment.experiment_id = signal_alignment.experiment_id 
                        and signal_alignment.is_complete=1 
                        and signal_alignment.has_duplicates=0
                        )
                        join data_file signal_bam on (signal_alignment.bam_file_id = signal_bam.data_file_id)
                        join experiment control_experiment on (control_experiment.experiment_id = experiment.control_id)
                        join alignment control_alignment on (
                        control_experiment.experiment_id = control_alignment.experiment_id 
                        and control_alignment.is_complete=1 
                        and control_alignment.has_duplicates=0
                        )
                        join data_file control_bam on (control_alignment.bam_file_id = control_bam.data_file_id)
                    where 
                        feature_type.name in (
                        "H3K4me1", 
                        "H3K4me2", 
                        "H3K4me3", 
                        "H3K9ac", 
                        "H3K9me3",
                        "H3K27ac", 
                        "H3K27me3", 
                        "H3K36me3", 
                        "DNase1",
                        "CTCF"
                        )
                    order by 
                        epigenome, feature_type
                ~,
                q~ drop table if exists temp_segmentation_epigenome_with_ctcf ~,
                q~
                    create table temp_segmentation_epigenome_with_ctcf
                    select 
                        epigenome
                    from 
                        temp_segmentation_candidates 
                    where 
                        feature_type in ("H3K4me1", "H3K4me3", "H3K27ac", "H3K27me3", "H3K36me3", "CTCF") 
                    group by 
                        epigenome 
                    having 
                        count(distinct feature_type) = 6
                ~,
                q~ drop table if exists temp_segmentation_epigenome_without_ctcf ~,
                q~
                    create table temp_segmentation_epigenome_without_ctcf
                    select 
                        distinct epigenome
                    from 
                        temp_segmentation_candidates
                        left join temp_segmentation_epigenome_with_ctcf using (epigenome)
                    where 
                        temp_segmentation_epigenome_with_ctcf.epigenome is null
                ~,
                q~ drop table if exists segmentation_cell_table_ctcf ~,
                q~
                    create table segmentation_cell_table_ctcf
                    select * from temp_segmentation_epigenome_with_ctcf join temp_segmentation_candidates using (epigenome)
                ~,
                q~ drop table if exists segmentation_cell_table_without_ctcf ~,
                q~
                    create table segmentation_cell_table_without_ctcf
                    select * from temp_segmentation_epigenome_without_ctcf join temp_segmentation_candidates using (epigenome)
                    limit 6
                ~,
                q~ drop table if exists temp_segmentation_candidates ~,
                q~ drop table if exists temp_segmentation_epigenome_with_ctcf ~,
                q~ drop table if exists temp_segmentation_epigenome_without_ctcf ~,
              ],
              db_conn => 'funcgen:#species#',
          },
          -flow_into => {
            MAIN => 'make_segmentation_dir',
          },
      },
      {   -logic_name => 'make_segmentation_dir',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
              cmd => 
                  qq( rm    -rf #tempdir_segmentation#/#species# ; )
                . qq( mkdir -p  #tempdir_segmentation#/#species#   ),
          },
          -flow_into => { 
            'MAIN->A' => 'start_parallel_segmentation',
            'A->MAIN' => 'start_regulatory_build',
          },
      },
      {   -logic_name => 'start_parallel_segmentation',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -flow_into => { 
            MAIN     => [
                'generate_cell_table_file_with_ctcf',
                'generate_cell_table_file_without_ctcf',
            ],
          },
      },
      {   -logic_name => 'generate_cell_table_file_with_ctcf',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => qq( generate_cell_table_file.pl      )
                 . qq(   --species         #species#    )
                 . qq(   --registry        #reg_conf#   )
                 . qq(   --cell_table      segmentation_cell_table_ctcf )
                 . qq(   --cell_table_file #tempdir_segmentation#/#species#/celltable.with_ctcf.txt )
          },
          -flow_into => {
            MAIN     => 'binarize_with_ctcf',
          },
      },
      {   -logic_name => 'generate_cell_table_file_without_ctcf',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => qq( generate_cell_table_file.pl      )
                 . qq(   --species         #species#    )
                 . qq(   --registry        #reg_conf#   )
                 . qq(   --cell_table      segmentation_cell_table_without_ctcf )
                 . qq(   --cell_table_file #tempdir_segmentation#/#species#/celltable.without_ctcf.txt )
          },
          -flow_into => {
            MAIN     => 'binarize_without_ctcf',
          },
      },
      {   -logic_name => 'binarize_with_ctcf',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -rc_name    => 'binarization',
          -parameters => { 
            cmd => qq( java -Xmx14500m                     )
                 . qq(   -jar #ChromHMM#                   )
                 . qq(   BinarizeBam                       )
                 . qq(   #chromosome_length_file#          )
                 . qq(   #data_root_dir_species_assembly#  )
                 . qq(   #tempdir_segmentation#/#species#/celltable.with_ctcf.txt )
                 . qq(   #binarized_bam_dir#/with_ctcf/               ),
          },
          -flow_into => {
            MAIN     => 'learn_model_with_ctcf',
          },
      },
      {   -logic_name => 'binarize_without_ctcf',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -rc_name    => 'binarization',
          -parameters => { 
            cmd => qq( java -Xmx14500m                     )
                 . qq(   -jar #ChromHMM#                   )
                 . qq(   BinarizeBam                       )
                 . qq(   #chromosome_length_file#          )
                 . qq(   #data_root_dir_species_assembly#  )
                 . qq(   #tempdir_segmentation#/#species#/celltable.without_ctcf.txt )
                 . qq(   #binarized_bam_dir#/without_ctcf/               ),
          },
          -flow_into => {
            MAIN     => 'learn_model_without_ctcf',
          },
      },
      {   -logic_name => 'learn_model_with_ctcf',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -rc_name    => 'learn_model',
          -parameters => {
            cmd => 
                qq( export DISPLAY= ;             )
                . qq( java -Xmx30000m             )
                . qq( -jar #ChromHMM#             )
                . qq( LearnModel                  )
                . qq( -r 50 -p 12                )
                . qq( -l #chromosome_length_file# )
                . qq( #binarized_bam_dir#/with_ctcf/         )
                . qq( #learnmodel_output_dir#/with_ctcf/     )
                . qq( 25                          )
                . qq( #assembly#                  )
          },
      },
      {   -logic_name => 'learn_model_without_ctcf',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -rc_name    => 'learn_model',
          -parameters => {
            cmd => 
                qq( export DISPLAY= ;             )
                . qq( java -Xmx30000m             )
                . qq( -jar #ChromHMM#             )
                . qq( LearnModel                  )
                . qq( -r 50 -p 12                )
                . qq( -l #chromosome_length_file# )
                . qq( #binarized_bam_dir#/without_ctcf/         )
                . qq( #learnmodel_output_dir#/without_ctcf/     )
                . qq( 25                          )
                . qq( #assembly#                  )
          },
      },
      
      {   -logic_name => 'start_regulatory_build',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -flow_into => { 
            MAIN => 'generate_segmentation_parameter_file',
          },
      },
      
      
      {   -logic_name => 'generate_segmentation_parameter_file',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => q( echo -e "segmentation\t#species#_segmentation_with_ctcf\tChromHMM\t#learnmodel_output_dir#/with_ctcf/"        > #segmentation_parameter_file#  ; )
                 . q( echo -e "segmentation\t#species#_segmentation_without_ctcf\tChromHMM\t#learnmodel_output_dir#/without_ctcf/" >> #segmentation_parameter_file#  ; )
          },
          -flow_into => { 
            MAIN     => 'make_regbuild_dir',
          },
      },
      {   -logic_name => 'make_regbuild_dir',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
              cmd => 
                  qq(rm    -rf #tempdir_regulatory_build#/#species#/#assembly# ; )
                . qq(mkdir -p  #tempdir_regulatory_build#/#species#/#assembly#   ),
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
            MAIN => [ 
                'load_build', 
                'load_state_emissions_and_assignments',
            ]
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
            MAIN => [
                'create_regulatory_build_statistics', 
                'populate_regulatory_build_epigenome_table',
#                 'register_segmentation_files',
            ]
          },
      },
      {   -logic_name => 'load_state_emissions_and_assignments',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 
                q(
                  load_state_emissions_and_assignments.pl \
                      --species          #species# \
                      --registry         #reg_conf# \
                      --emissions_file   #learnmodel_output_dir#/with_ctcf/emissions_25.txt \
                      --segmentation     #species#_segmentation_with_ctcf \
                      --assignments_file #tempdir_regulatory_build#/#species#/tmp/assignments.txt
                )
          },
          -flow_into => {
            MAIN => 'load_state_emissions_without_ctcf',
          },
      },
      {   -logic_name => 'load_state_emissions_without_ctcf',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 
                q(
                  load_state_emissions_and_assignments.pl \
                      --species          #species# \
                      --registry         #reg_conf# \
                      --emissions_file   #learnmodel_output_dir#/without_ctcf/emissions_25.txt \
                      --segmentation     #species#_segmentation_without_ctcf
                )
          },
          -flow_into => {
            MAIN => 'generate_segmentation_report',
          },
      },
      {   -logic_name => 'generate_segmentation_report',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 
                q(
                  generate_segmentation_report.pl \
                    --species          #species# \
                    --registry         #reg_conf# \
                    --output_directory #reports_dir#/#species#
                )
          },
#           -flow_into => {
#             MAIN => 'create_regulatory_build_statistics',
#           },
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
                "
                insert into regulatory_build_statistics (regulatory_build_id, statistic, value) (
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
                insert into regulatory_build_statistics (regulatory_build_id, statistic, value) (
                    select regulatory_build_id, 'number_promoter', count(regulatory_feature_id) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "Promoter" group by feature_type.name, regulatory_build_id
                );
                ~,
                qq~
                insert into regulatory_build_statistics (regulatory_build_id, statistic, value) (
                    select regulatory_build_id, 'number_enhancer', count(regulatory_feature_id) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "Enhancer" group by feature_type.name, regulatory_build_id
                );
                ~,
                qq~
                insert into regulatory_build_statistics (regulatory_build_id, statistic, value) (
                    select regulatory_build_id, 'number_promoter_flanking_region', count(regulatory_feature_id) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "Promoter Flanking Region" group by feature_type.name, regulatory_build_id
                );
                ~,
                qq~
                insert into regulatory_build_statistics (regulatory_build_id, statistic, value) (
                    select regulatory_build_id, 'number_transcription_factor_binding_site', count(regulatory_feature_id) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "TF binding site" group by feature_type.name, regulatory_build_id
                );
                ~,
                qq~
                insert into regulatory_build_statistics (regulatory_build_id, statistic, value) (
                    select regulatory_build_id, 'number_open_chromatin', count(regulatory_feature_id) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "Open chromatin" group by feature_type.name, regulatory_build_id
                );
                ~,
                qq~
                insert into regulatory_build_statistics (regulatory_build_id, statistic, value) (
                    select regulatory_build_id, 'number_ctcf_binding_site', count(regulatory_feature_id) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "CTCF Binding Site" group by feature_type.name, regulatory_build_id
                );
                ~,

                qq~
                insert into regulatory_build_statistics (regulatory_build_id, statistic, value) (
                    select regulatory_build_id, 'average_length_promoter', AVG( (seq_region_end + bound_end_length ) - (seq_region_start-bound_start_length)  + 1) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "Promoter" group by feature_type.name, regulatory_build_id
                );
                ~,
                qq~
                insert into regulatory_build_statistics (regulatory_build_id, statistic, value) (
                    select regulatory_build_id, 'average_length_enhancer', AVG( (seq_region_end + bound_end_length ) - (seq_region_start-bound_start_length)  + 1) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "Enhancer" group by feature_type.name, regulatory_build_id
                );
                ~,
                qq~
                insert into regulatory_build_statistics (regulatory_build_id, statistic, value) (
                    select regulatory_build_id, 'average_length_promoter_flanking_region', AVG( (seq_region_end + bound_end_length ) - (seq_region_start-bound_start_length)  + 1) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "Promoter Flanking Region" group by feature_type.name, regulatory_build_id
                );
                ~,
                qq~
                insert into regulatory_build_statistics (regulatory_build_id, statistic, value) (
                    select regulatory_build_id, 'average_length_transcription_factor_binding_site', AVG( (seq_region_end + bound_end_length ) - (seq_region_start-bound_start_length)  + 1) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "TF binding site" group by feature_type.name, regulatory_build_id
                );
                ~,
                qq~
                insert into regulatory_build_statistics (regulatory_build_id, statistic, value) (
                    select regulatory_build_id, 'average_length_open_chromatin', AVG( (seq_region_end + bound_end_length ) - (seq_region_start-bound_start_length)  + 1) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "Open chromatin" group by feature_type.name, regulatory_build_id
                );
                ~,
                qq~
                insert into regulatory_build_statistics (regulatory_build_id, statistic, value) (
                    select regulatory_build_id, 'average_length_ctcf_binding_site', AVG( (seq_region_end + bound_end_length ) - (seq_region_start-bound_start_length)  + 1) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "CTCF Binding Site" group by feature_type.name, regulatory_build_id
                );
                ~,

                qq~
                insert into regulatory_build_statistics (regulatory_build_id, statistic, value) (
                    select regulatory_build_id, 'sum_length_promoter', SUM( (seq_region_end + bound_end_length ) - (seq_region_start-bound_start_length)  + 1) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "Promoter" group by feature_type.name, regulatory_build_id
                );
                ~,
                qq~
                insert into regulatory_build_statistics (regulatory_build_id, statistic, value) (
                    select regulatory_build_id, 'sum_length_enhancer', SUM( (seq_region_end + bound_end_length ) - (seq_region_start-bound_start_length)  + 1) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "Enhancer" group by feature_type.name, regulatory_build_id
                );
                ~,
                qq~
                insert into regulatory_build_statistics (regulatory_build_id, statistic, value) (
                    select regulatory_build_id, 'sum_length_promoter_flanking_region', SUM( (seq_region_end + bound_end_length ) - (seq_region_start-bound_start_length)  + 1) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "Promoter Flanking Region" group by feature_type.name, regulatory_build_id
                );
                ~,
                qq~
                insert into regulatory_build_statistics (regulatory_build_id, statistic, value) (
                    select regulatory_build_id, 'sum_length_transcription_factor_binding_site', SUM( (seq_region_end + bound_end_length ) - (seq_region_start-bound_start_length)  + 1) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "TF binding site" group by feature_type.name, regulatory_build_id
                );
                ~,
                qq~
                insert into regulatory_build_statistics (regulatory_build_id, statistic, value) (
                    select regulatory_build_id, 'sum_length_open_chromatin', SUM( (seq_region_end + bound_end_length ) - (seq_region_start-bound_start_length)  + 1) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "Open chromatin" group by feature_type.name, regulatory_build_id
                );
                ~,
                qq~
                insert into regulatory_build_statistics (regulatory_build_id, statistic, value) (
                    select regulatory_build_id, 'sum_length_ctcf_binding_site', SUM( (seq_region_end + bound_end_length ) - (seq_region_start-bound_start_length)  + 1) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "CTCF Binding Site" group by feature_type.name, regulatory_build_id
                );
                ~

              ],
              db_conn => 'funcgen:#species#',
          },
          -flow_into => {
            MAIN => 'generate_regulatory_build_report',
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
            MAIN => 'register_segmentation_files',
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
                      --projected_segmentation_directory #tempdir_regulatory_build#/#species#/#assembly#/projected_segmentations \
                      --db_file_species_assembly_dir     #data_root_dir#/#species#/#assembly#/funcgen/segmentation_file/#ensembl_release_version# \
                      --db_file_relative_dir             /funcgen/segmentation_file/#ensembl_release_version#
                )
          },
#           -flow_into => {
#             MAIN => '',
#           },
      },
      {   -logic_name => 'generate_regulatory_build_report',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 
                q(
                  generate_regulatory_build_report.pl  \
                      --species          #species# \
                      --registry         #reg_conf# \
                      --output_directory #reports_dir#/#species#
                )
          },
#           -flow_into => {
#             MAIN => 'segmentation_done',
#           },
      },
#       {   -logic_name => 'segmentation_done',
#           -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
#       },
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
