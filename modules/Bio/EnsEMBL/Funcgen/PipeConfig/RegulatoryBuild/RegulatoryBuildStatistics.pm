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
                'truncate regulatory_build_statistic;',
              ],
              db_conn => 'funcgen:#species#',
          },
          -flow_into   => {
              MAIN => 'compute_genome_coverage',
          },
      },
      {   -logic_name => 'compute_genome_coverage',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 
                q(
                  compute_genome_coverage.pl \
                    --species  #species# \
                    --registry #reg_conf# \
                    --genome_coverage_report #tempdir_regulatory_build_statistics#/#species#/genome_coverage_report.pl
                )
          },
          -flow_into => {
            MAIN => 'store_genome_coverage',
          },
      },
      {   -logic_name => 'store_genome_coverage',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 
                q(
                    store_genome_coverage_report.pl \
                        --species  #species# \
                        --registry #reg_conf# \
                        --genome_coverage_report #tempdir_regulatory_build_statistics#/#species#/genome_coverage_report.pl
                )
          },
          -flow_into => {
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
                    --tempdir  #tempdir_regulatory_build_statistics#/#species#
                )
          },

          -flow_into => {
            MAIN => 'compute_rb_basic_statistics',
          },
      },
      {   -logic_name => 'compute_rb_basic_statistics',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 
                q(
                  compute_regulatory_build_basic_statistics.pl \
                    --species  #species# \
                    --registry #reg_conf# \
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
                  compute_regulatory_build_quantiles.pl \
                    --species  #species# \
                    --registry #reg_conf# \
                )
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
                      --output_file      #reports_dir#/#species#/regulatory_build.html
                )
          },
      },
    ];
}

1;
