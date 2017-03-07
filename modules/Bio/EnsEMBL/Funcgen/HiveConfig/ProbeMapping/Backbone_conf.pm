package Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::Backbone_conf;

use strict;
use warnings;
use base ('Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::Base');
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub pipeline_analyses {
    my $self = shift;

    return [
        {   -logic_name  => 'start',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               MAIN => 'backbone_fire_import',
            },
        },
        {   -logic_name  => 'backbone_fire_import',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_import',
               'A->1' => 'backbone_fire_export'
            },
        },
        {
            -logic_name  => 'start_import',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        {   -logic_name  => 'backbone_fire_export',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_export',
               'A->1' => 'backbone_fire_align_probes'
            },
        },
        {
            -logic_name  => 'start_export',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        {   -logic_name  => 'backbone_fire_align_probes',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_align_probes',
               'A->1' => 'backbone_fire_probe2transcript'
            },
        },
        {
            -logic_name  => 'start_align_probes',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        {   -logic_name  => 'backbone_fire_probe2transcript',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_probe2transcript',
               'A->1' => 'backbone_fire_healthchecks'
            },
        },
        {
            -logic_name  => 'start_probe2transcript',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        {   -logic_name  => 'backbone_fire_healthchecks',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_healthchecks',
               'A->1' => 'backbone_pipeline_finished'
            },
        },
        {
            -logic_name  => 'start_healthchecks',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        {   -logic_name => 'backbone_pipeline_finished',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        }
    ]
}

1;
