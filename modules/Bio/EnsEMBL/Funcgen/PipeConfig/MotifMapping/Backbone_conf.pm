=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::PipeConfig::MotifMapping::Backbone_conf

=head1 DESCRIPTION



=head1 LICENSE

    Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
    Copyright [2016-2024] EMBL-European Bioinformatics Institute

    Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

         http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software distributed under the License
    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and limitations under the License.

=head1 CONTACT

    ensembl-funcgen@ebi.ac.uk

=cut

package Bio::EnsEMBL::Funcgen::PipeConfig::MotifMapping::Backbone_conf;

use strict;
use warnings;
use base ('Bio::EnsEMBL::Funcgen::PipeConfig::MotifMapping::Base_conf');
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub pipeline_analyses {
    my $self = shift;

    return [
        {
            -logic_name => 'start',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into  => {
                MAIN => [ 'backbone_fire_pre_pipeline_checks_tasks' ],
            },
        },
        {
            -logic_name => 'backbone_fire_pre_pipeline_checks_tasks',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into  => {
                '1->A' => 'start_pre_pipeline_checks_tasks',
                'A->1' => 'backbone_fire_start_registration',
            },
        },
        {
            -logic_name => 'start_pre_pipeline_checks_tasks',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        {
            -logic_name => 'backbone_fire_start_registration',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into  => {
                '1->A' => 'start_registration',
                'A->1' => 'backbone_junction'
            },
        },
        {
            -logic_name => 'start_registration',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        {
            -logic_name => 'backbone_junction',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into  => {
                '1->A' => [ 'backbone_fire_start_MOODS',
                            'backbone_fire_start_Peak_RF_BED_Dump' ],
                'A->1' => 'backbone_fire_start_PeaksIntersection'
            },
        },
        {
            -logic_name => 'backbone_fire_start_Peak_RF_BED_Dump',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into  => {
                '1' => 'generate_Peak_RF_BED_Files',
            },
        },
        {
            -logic_name => 'generate_Peak_RF_BED_Files',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },

        {
            -logic_name => 'backbone_fire_start_MOODS',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into  => {
                '1' => 'start_MOODS',
            },
        },
        {
            -logic_name => 'start_MOODS',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        {
            -logic_name => 'backbone_fire_start_PeaksIntersection',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into  => {
                '1->A' => 'start_PeaksIntersection',
                'A->1' => 'backbone_fire_start_RegulatoryFeatureIntersection'
            }
        },
        {
            -logic_name => 'start_PeaksIntersection',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',

        },
        {
            -logic_name => 'backbone_fire_start_RegulatoryFeatureIntersection',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into  => {
                '1->A' => 'start_RegulatoryFeaturesIntersection',
                'A->1' => 'backbone_fire_start_StableIdMapping'
            }
        },
        {
            -logic_name => 'start_RegulatoryFeaturesIntersection',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',

        },
        {
            -logic_name => 'backbone_fire_start_StableIdMapping',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into  => {
                '1->A' => 'start_StableIdMapping',
                'A->1' => 'backbone_fire_start_WebBedFile'
            }
        },
        {
            -logic_name => 'start_StableIdMapping',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',

        },
        {
            -logic_name => 'backbone_fire_start_WebBedFile',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into  => {
                '1->A' => 'start_WebBedFile',
                'A->1' => 'backbone_fire_start_FtpDumps'
            }
        },
        {
            -logic_name => 'start_WebBedFile',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',

        },
        {
            -logic_name => 'backbone_fire_start_FtpDumps',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into  => {
                '1->A' => 'start_FtpDumps',
                'A->1' => 'end'
            }
        },
        {
            -logic_name => 'start_FtpDumps',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',

        },
        {
            -logic_name => 'end',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        }
    ];
}

1;