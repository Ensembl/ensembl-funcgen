package Bio::EnsEMBL::Funcgen::PipeConfig::ChIPSeqAnalysis::Idr;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Funcgen::PipeConfig::ChIPSeqAnalysis::Base';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub pipeline_analyses {
    my $self = shift;
    return [
        {   -logic_name  => 'start_idr',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               MAIN => [
                 WHEN('#execution_plan#->{idr}->{strategy} ne "SKIP_IDR"' => ['start_the_idr']),
                 WHEN('#execution_plan#->{idr}->{strategy} eq "SKIP_IDR"' => ['store_no_idr'])
               ]
            },
        },
       {   -logic_name  => 'store_no_idr',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::StoreNoIdr',
        },
       {   -logic_name  => 'start_the_idr',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               'MAIN->A' => 'seed_swembl_runs',
               'A->MAIN' => 'compute_idr_per_experiment',
            },
        },
        {   -logic_name  => 'compute_idr_per_experiment',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::ComputeIdrPerExperiment',
            -flow_into   => {
               '2->A' => 'seed_pairwise',
               'A->2' => 'idr_peak_threshold',
            },
        },
        {   -logic_name  => 'seed_swembl_runs',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::SeedSwemblRuns',
            -flow_into   => {
               2 => 'run_swembl',
            },
        },
        {   -logic_name  => 'run_swembl',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::RunPermissiveSWEmbl',
            -flow_into   => {
               2 => [
                #'?accu_name=permissive_peak_file&accu_address=[]&accu_input_variable=permissive_peak_file',
                '?accu_name=permissive_peak_calling&accu_address=[]&accu_input_variable=permissive_peak_calling',
               ],
            },
        },
        {   -logic_name  => 'seed_pairwise',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::SeedPairwise',
            -flow_into   => {
               2 => 'run_idr'
            },
        },
        {   -logic_name  => 'run_idr',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::RunIDR',
            -flow_into   => {
               2 => '?accu_name=idr_result&accu_address=[]&accu_input_variable=idr_result',
            }
        },
        {   -logic_name  => 'idr_peak_threshold',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::IdrPeakThreshold',
        },
    ];
}

1;
