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
               MAIN => 'pre_pipeline_checks',
            },
        },
        {   -logic_name  => 'pre_pipeline_checks',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::PrePipelineChecks',
            -flow_into => {
                MAIN => 'make_temp_dir',
            },
        },
        {
            -logic_name  => 'make_temp_dir',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd       => 'mkdir -p #tempdir#/#species#',
            },
            -flow_into => {
                MAIN => 'backbone_fire_import',
            },
        },
        {   -logic_name  => 'backbone_fire_import',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_import',
               'A->1' => 'backbone_fire_import_healthchecks'
            },
        },
        {
            -logic_name  => 'start_import',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        {   -logic_name  => 'backbone_fire_import_healthchecks',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_import_healthchecks',
               'A->1' => 'backbone_fire_export'
            },
        },
        {
            -logic_name  => 'start_import_healthchecks',
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
               'A->1' => 'backbone_fire_align_healthchecks'
            },
        },
        {
            -logic_name  => 'start_align_probes',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        {   -logic_name  => 'backbone_fire_align_healthchecks',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_align_healthchecks',
               'A->1' => 'backbone_fire_probe2transcript'
            },
        },
        {
            -logic_name  => 'start_align_healthchecks',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        {   -logic_name  => 'backbone_fire_probe2transcript',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_probe2transcript',
               'A->1' => 'backbone_fire_probe_to_transcript_healthchecks'
            },
        },
        {
            -logic_name  => 'start_probe2transcript',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },

        {   -logic_name  => 'backbone_fire_probe_to_transcript_healthchecks',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_probe_to_transcript_healthchecks',
               'A->1' => 'backbone_fire_switch_table_engines'
            },
        },
        {
            -logic_name  => 'start_probe_to_transcript_healthchecks',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        {   -logic_name  => 'backbone_fire_switch_table_engines',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_switch_table_engines',
               'A->1' => 'backbone_fire_switch_table_engines_healthchecks'
            },
        },
        {
            -logic_name  => 'start_switch_table_engines',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },

        {   -logic_name  => 'backbone_fire_switch_table_engines_healthchecks',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_switch_table_engine_healthchecks',
               'A->1' => 'backbone_pipeline_finished'
            },
        },
        {
            -logic_name  => 'start_switch_table_engine_healthchecks',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },

        {   -logic_name => 'backbone_pipeline_finished',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        }


    ]
}

1;
