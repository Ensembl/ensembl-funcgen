package Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::PrePipelineChecks;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::Base';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub pipeline_analyses {
    my $self = shift;

    return [
        {   -logic_name  => 'start_pre_pipeline_checks',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               'MAIN->A' => 'pre_pipeline_checks',
               'A->MAIN' => 'truncate_execution_plan_table',
            },
        },
        {
          -logic_name => 'pre_pipeline_checks',
          -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::PrePipelineChecks',
          -max_retry_count => 0,
        },
        {
            -logic_name  => 'truncate_execution_plan_table',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
            -parameters => {
                sql     => [
                    "truncate execution_plan;",
                ],
                db_conn => 'funcgen:#species#',
            },
            -flow_into   => {
               'MAIN->A' => 'fetch_experiments_to_process',
               'A->MAIN' => 'check_execution_plans',
            },
        },
        {   -logic_name  => 'fetch_experiments_to_process',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -parameters => {
                db_conn    => 'funcgen:#species#',
                
                # Everything
                inputquery => '
                  select 
                    experiment_id,
                    experiment.name as experiment_name
                  from 
                    experiment 
                    join feature_type using (feature_type_id) 
                    join experimental_group using (experimental_group_id) 
                  where 
                    is_control = 0 
                    and class in (
                      "Histone", 
                      "Open Chromatin",
                      "Transcription Factor",
                      "Polymerase"
                    )
                ',
#                  and (experiment.name like "K562%")
#                    and (experiment.name like "Lung%" or experiment.name like "HepG2%")

            },
            -flow_into => {
               2 => { 'create_execution_plan' => INPUT_PLUS },
            },
        },
        {   -logic_name  => 'create_execution_plan',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::CreateExecutionPlan',
            -analysis_capacity => 20,
        },
        {   -logic_name  => 'check_execution_plans',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::CheckExecutionPlans',
            -flow_into   => {
               1 => 'backbone_fire_start_check_paired_end_read_files'
            },
        },
        {   -logic_name  => 'backbone_fire_start_check_paired_end_read_files',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               MAIN => 'start_check_paired_end_read_files',
            },
        },
        {   -logic_name => 'start_check_paired_end_read_files',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
    ]
}

1;
