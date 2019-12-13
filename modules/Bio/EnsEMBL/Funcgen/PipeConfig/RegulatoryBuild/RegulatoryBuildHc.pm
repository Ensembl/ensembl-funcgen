=pod 

=head1 NAME

    Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::ChIPSeqCleanup

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

package Bio::EnsEMBL::Funcgen::PipeConfig::RegulatoryBuild::RegulatoryBuildHc;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base 'Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::Base';

sub pipeline_analyses {
    my ($self) = @_;
    return [
        {   -logic_name => 'start_regulatory_build_hc',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into => { 
              MAIN => [
                'hc_current_regulatory_build_exists',
                'hc_regulatory_build_epigenomes_populated',
                'hc_regulatory_activities_counts',
                'hc_segmentation_files_exist',
                'regulatory_build_checks',
              ]
            },
        },
        {
          -logic_name       => 'regulatory_build_checks',
          -module           => 'Bio::EnsEMBL::DataCheck::Pipeline::RunDataChecks',
          -max_retry_count  => 0,
          -parameters => {
            registry_file    => '#reg_conf#',
            species          => '#species#',
            group            => 'funcgen',
            datacheck_groups => [ 'regulatory_build' ],
          },
        },
        {   -logic_name => 'hc_segmentation_files_exist',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
              cmd => 
                  q(
                    check_segmentation_files_exist.pl \
                      --species      #species#        \
                      --registry     #reg_conf#       \
                      --db_file_path #data_root_dir#
                  )
            },
        },
        {
            -logic_name  => 'hc_current_regulatory_build_exists',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SqlHealthcheck',
            -parameters => {
              db_conn       => 'funcgen:#species#',
              description   => 'Make sure there is a current regulatory build',
              query         => "
                select regulatory_build_id from regulatory_build where is_current = true
              ",
              expected_size => '1'
            },
        },
        {
            -logic_name  => 'hc_regulatory_build_epigenomes_populated',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SqlHealthcheck',
            -parameters => {
              db_conn       => 'funcgen:#species#',
              description   => 'Make sure the epigenomes present in the regulatory build have been registered in the dedicated table',
              query         => "
                select 
                  regulatory_build_id, 
                  count(epigenome_id) as num_epigenomes 
                from 
                  regulatory_build 
                  left join regulatory_build_epigenome using (regulatory_build_id) 
                having 
                  num_epigenomes = 0
              ",
              expected_size => '0'
            },
        },
        {
            -logic_name  => 'hc_regulatory_activities_counts',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SqlHealthcheck',
            -parameters => {
              db_conn       => 'funcgen:#species#',
              description   => 'All regulatory features should have as many regulatory activities as there are epigenomes in the regulatory build.',
              query         => "
                select 
                  regulatory_feature.regulatory_feature_id, 
                  count(regulatory_activity_id) as num_activities 
                from 
                  regulatory_build 
                  join regulatory_feature on (
                    regulatory_build.regulatory_build_id = regulatory_feature.regulatory_build_id 
                    and regulatory_build.is_current = true
                  ) 
                  left join regulatory_activity using (regulatory_feature_id) 
                group by 
                  regulatory_feature_id 
                having 
                  num_activities not in (
                    select 
                      count(epigenome_id) 
                    from 
                      regulatory_build_epigenome 
                      join regulatory_build 
                    where 
                      is_current = true
                );
              ",
              expected_size => '0'
            },
        },
    ];
}

1;
