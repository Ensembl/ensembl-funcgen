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
package Bio::EnsEMBL::Funcgen::PipeConfig::RegulatoryBuild::RegulatoryBuildStatistics;

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
    tempdir_regulatory_build_statistics => $self->o('tempdir') . '/regulatory_build_statistics',
    data_root_dir            => $self->o('data_root_dir'),
    reference_data_root_dir  => $self->o('reference_data_root_dir'),
    reg_conf                 => $self->o('reg_conf'),
  };
}

sub pipeline_analyses {
    my ($self) = @_;
    return [
      {   -logic_name => 'start_regulatory_build_statistics',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -flow_into => { 
            MAIN => 'truncate_regulatory_build_statistics_tables',
          },
      },
      {
          -logic_name  => 'truncate_regulatory_build_statistics_tables',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
          -parameters => {
              sql     => [
                "truncate regulatory_build_statistic;",
              ],
              db_conn => 'funcgen:#species#',
          },
          -flow_into   => {
              MAIN => 'compute_enhancer_coverage',
          },
      },
      {   -logic_name => 'compute_enhancer_coverage',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
#             cmd => 
#                 q(
#                   compute_enhancer_coverage.pl \
#                     --species  #species# \
#                     --registry #reg_conf#
#                 )
#           },
            cmd => 
                q(
                  compute_enhancer_coverage_bedtools.pl \
                    --species  #species# \
                    --registry #reg_conf# \
                    --tempdir #tempdir_regulatory_build_statistics#/#species#
                )
          },

          -flow_into => {
            MAIN => 'compute_quantiles',
          },
      },
      {   -logic_name => 'compute_quantiles',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 
                q(
                  compute_regulatory_build_statistics.pl \
                    --species  #species# \
                    --registry #reg_conf# \
                )
          },

          -flow_into => {
            MAIN => 'create_regulatory_build_statistic',
          },
      },
      {
          -logic_name  => 'create_regulatory_build_statistic',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
          -parameters => {
              sql     => [
                "drop table if exists regulatory_build_statistic",

                "CREATE TABLE regulatory_build_statistic (
                  regulatory_build_statistic_id INT(30) UNSIGNED NOT NULL AUTO_INCREMENT,
                  regulatory_build_id INT(22) UNSIGNED NOT NULL,
                  statistic                VARCHAR(128) NOT NULL,
                  value                    BIGINT(11) UNSIGNED DEFAULT '0' NOT NULL,
                  PRIMARY KEY (regulatory_build_statistic_id),
                  UNIQUE KEY stats_uniq(statistic, regulatory_build_id)
                )",
                "insert into regulatory_build_statistic (regulatory_build_id, statistic, value) (
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

#                 qq~
#                 insert into regulatory_build_statistic (regulatory_build_id, statistic, value) (
#                     select regulatory_build_id, 'average_length_promoter', AVG( (seq_region_end + bound_end_length ) - (seq_region_start-bound_start_length)  + 1) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "Promoter" group by feature_type.name, regulatory_build_id
#                 );
#                 ~,
#                 qq~
#                 insert into regulatory_build_statistic (regulatory_build_id, statistic, value) (
#                     select regulatory_build_id, 'average_length_enhancer', AVG( (seq_region_end + bound_end_length ) - (seq_region_start-bound_start_length)  + 1) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "Enhancer" group by feature_type.name, regulatory_build_id
#                 );
#                 ~,
#                 qq~
#                 insert into regulatory_build_statistic (regulatory_build_id, statistic, value) (
#                     select regulatory_build_id, 'average_length_promoter_flanking_region', AVG( (seq_region_end + bound_end_length ) - (seq_region_start-bound_start_length)  + 1) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "Promoter Flanking Region" group by feature_type.name, regulatory_build_id
#                 );
#                 ~,
#                 qq~
#                 insert into regulatory_build_statistic (regulatory_build_id, statistic, value) (
#                     select regulatory_build_id, 'average_length_transcription_factor_binding_site', AVG( (seq_region_end + bound_end_length ) - (seq_region_start-bound_start_length)  + 1) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "TF binding site" group by feature_type.name, regulatory_build_id
#                 );
#                 ~,
#                 qq~
#                 insert into regulatory_build_statistic (regulatory_build_id, statistic, value) (
#                     select regulatory_build_id, 'average_length_open_chromatin', AVG( (seq_region_end + bound_end_length ) - (seq_region_start-bound_start_length)  + 1) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "Open chromatin" group by feature_type.name, regulatory_build_id
#                 );
#                 ~,
#                 qq~
#                 insert into regulatory_build_statistic (regulatory_build_id, statistic, value) (
#                     select regulatory_build_id, 'average_length_ctcf_binding_site', AVG( (seq_region_end + bound_end_length ) - (seq_region_start-bound_start_length)  + 1) from regulatory_feature join feature_type using (feature_type_id) where feature_type.name = "CTCF Binding Site" group by feature_type.name, regulatory_build_id
#                 );
#                 ~,

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

              ],
              db_conn => 'funcgen:#species#',
          },
          -flow_into => {
            MAIN => 'generate_regulatory_build_report',
          },
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
      },
    ];
}

1;
