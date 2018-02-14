package Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::Idr;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::Base';
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
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::StoreNoIdr',
        },
       {   -logic_name  => 'start_the_idr',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               'MAIN->A' => 'seed_swembl_runs',
               'A->MAIN' => 'compute_idr_per_experiment',
            },
        },
        {   -logic_name  => 'compute_idr_per_experiment',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::ComputeIdrPerExperiment',
            -flow_into   => {
               '2->A' => 'seed_pairwise',
               'A->2' => 'idr_peak_threshold',
            },
        },
        {   -logic_name  => 'seed_swembl_runs',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::SeedSwemblRuns',
            -flow_into   => {
               2 => 'run_swembl',
            },
        },
        {   -logic_name  => 'run_swembl',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::RunPermissiveSWEmbl',
            -parameters => {
              tempdir => '#tempdir_peak_calling#/#species#/idr',
            },
            -flow_into   => {
               2 => '?accu_name=permissive_peak_calling&accu_address=[]&accu_input_variable=permissive_peak_calling',
               MEMLIMIT => 'run_swembl_himem',
            },
        },
        {   -logic_name  => 'run_swembl_himem',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::RunPermissiveSWEmbl',
            -parameters => {
              tempdir => '#tempdir_peak_calling#/#species#/idr',
            },
            -rc_name     => '8Gb_job',
            -flow_into   => {
               2 => '?accu_name=permissive_peak_calling&accu_address=[]&accu_input_variable=permissive_peak_calling',
            },
        },
        {   -logic_name  => 'seed_pairwise',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::SeedPairwise',
            -flow_into   => {
               2 => 'run_idr'
            },
        },
        {   -logic_name  => 'run_idr',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::RunIDR',
            -flow_into   => {
               2        => '?accu_name=idr_result&accu_address=[]&accu_input_variable=idr_result',
               MEMLIMIT => 'run_idr_himem'
            }
        },
        {   -logic_name  => 'run_idr_himem',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::RunIDR',
            -rc_name => '8Gb_job',
            -flow_into   => {
               2 => '?accu_name=idr_result&accu_address=[]&accu_input_variable=idr_result',
            }
        },
        {   -logic_name  => 'idr_peak_threshold',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::IdrPeakThreshold',
        },
    ];
}

1;

