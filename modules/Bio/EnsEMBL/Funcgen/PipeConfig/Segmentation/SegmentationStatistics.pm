=pod 

=head1 NAME

    Bio::EnsEMBL::Funcgen::PipeConfig::Segmentation

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 LICENSE

    Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
    Copyright [2016-2020] EMBL-European Bioinformatics Institute

    Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

         http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software distributed under the License
    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and limitations under the License.

=head1 CONTACT

=cut
package Bio::EnsEMBL::Funcgen::PipeConfig::Segmentation::SegmentationStatistics;

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
    data_root_dir            => $self->o('data_root_dir'),
    reg_conf                 => $self->o('reg_conf'),
  };
}

sub pipeline_analyses {
    my ($self) = @_;
    return [
      {   -logic_name => 'start_segmentation_statistics',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -flow_into => { 
            MAIN => 'truncate_segmentation_statistic',
          },
      },
      {
          -logic_name  => 'truncate_segmentation_statistic',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
          -parameters => {
              sql     => [
                'truncate segmentation_statistic'
              ],
              db_conn => 'funcgen:#species#',
          },
          -flow_into   => {
              MAIN => 'load_segmentation_files_to_db',
          },
      },
      {   -logic_name => 'load_segmentation_files_to_db',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 
                q(
                  load_segmentation_files_to_db.pl     \
                    --species          #species#       \
                    --registry         #reg_conf#      \
                    --db_file_path     #data_root_dir# \
                    --tempdir          #tempdir_segmentation#/#species#/segmentation_file_loading
                )
          },
          -flow_into   => {
            MAIN  => 'compute_basic_segmentation_statistics',
          },
      },
      {   -logic_name => 'compute_basic_segmentation_statistics',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 
                q(
                    compute_basic_segmentation_statistics.pl \
                      --species  #species#             \
                      --registry #reg_conf#
                )
          },
          -flow_into   => {
            MAIN  => 'compute_segmentation_statistics',
          },
      },
      {   -logic_name => 'compute_segmentation_statistics',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -rc_name    => '32Gb_job',
          -parameters => {
            cmd => 
                q(
                    compute_segmentation_statistics.pl \
                      --species  #species#             \
                      --registry #reg_conf#
                )
          },
          -flow_into   => {
            MAIN  => 'generate_segmentation_report',
          },
      },
      {   -logic_name => 'generate_segmentation_report',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 
                q(
                  generate_segmentation_report.pl \
                    --species     #species# \
                    --registry    #reg_conf# \
                    --output_file #reports_dir#/#species#/segmentation.html
                )
          },
      },
    ];
}

1;
