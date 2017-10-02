package Bio::EnsEMBL::Funcgen::PipeConfig::ChIPSeqAnalysis;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Funcgen::PipeConfig::ChIPSeqAnalysis::Base';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub pipeline_analyses {
    my $self = shift;

    return [
        {   -logic_name  => 'start_chip_seq_analysis',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               'MAIN->A' => 'fetch_experiment_ids',
               'A->MAIN' => 'find_control_experiments',
            },
        },
        {   -logic_name  => 'fetch_experiment_ids',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -parameters => {
                db_conn    => 'funcgen:#species#',
                inputquery => '
                  select 
                    experiment_id 
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
                    and experimental_group.name != "BLUEPRINT"
                ',
            },
            -flow_into => {
               2 => { 'create_execution_plan' => INPUT_PLUS },
            },
        },
        {   -logic_name  => 'create_execution_plan',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::CreateExecutionPlan',
            -flow_into   => {
               2 => '?accu_name=execution_plan_list&accu_address=[]&accu_input_variable=plan',
            },
        },
        {   -logic_name  => 'find_control_experiments',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::FindControlExperiments',
            -flow_into   => {
               'A->2' => 'seed_bigwig_control_jobs',
               '3->A' => 'start_align_controls',
            },
        },
        {   -logic_name  => 'start_align_controls',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        {   -logic_name  => 'seed_bigwig_control_jobs',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               MAIN => [ 
                'seed_signal_processing',
                'write_bigwig'
              ],
            },
        },
        {   -logic_name  => 'write_bigwig',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::WriteBigWig',
        },
        {   -logic_name  => 'seed_signal_processing',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::SeedSignalProcessing',
            -flow_into   => {
               '2->A' => 'split_by_idr_strategy',
               'A->3' => 'call_peaks',
            },
        },
        {   -logic_name  => 'call_peaks',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               MAIN => 'done_chip_seq_analysis',
            },
        },
        {   -logic_name  => 'split_by_idr_strategy',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::SplitByIdrStrategy',
            -flow_into   => {
               '2->A' => 'no_idr',
               'A->2' => 'no_idr_done',
               '3->B' => 'idr_on_technical_replicates',
               'B->3' => 'idr_on_technical_replicates_done',
               '4->C' => 'idr_on_biological_replicates',
               'C->4' => 'idr_on_biological_replicates_done',
            },
        },
        {   -logic_name => 'done_chip_seq_analysis',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },

        #
        # idr_on_biological_replicates
        #
        {   -logic_name  => 'idr_on_biological_replicates',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        {   -logic_name => 'idr_on_biological_replicates_done',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        
        #
        # idr_on_technical_replicates
        #
        {   -logic_name  => 'idr_on_technical_replicates',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        {   -logic_name => 'idr_on_technical_replicates_done',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },

        #
        # no_idr
        #
        {   -logic_name  => 'no_idr',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        {   -logic_name => 'no_idr_done',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
    ]
}

1;

