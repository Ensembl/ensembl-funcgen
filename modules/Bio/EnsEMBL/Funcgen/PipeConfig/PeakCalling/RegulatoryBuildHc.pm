=pod 

=head1 NAME

    Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::ChIPSeqCleanup

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

package Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::RegulatoryBuildHc;

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
              MAIN => 'generate_peak_calling_report',
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
                      is_current=true
                );
              ",
              expected_size => '0'
            },
          -flow_into => {
              MAIN => 'regulatory_build_hc_done',
          },
        },
        {   -logic_name => 'regulatory_build_hc_done',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
    ];
}

1;
