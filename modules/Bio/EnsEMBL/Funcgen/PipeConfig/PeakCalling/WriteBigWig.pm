package Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::WriteBigWig;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::Base';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub pipeline_analyses {
    my $self = shift;
    return [
        {   -logic_name  => 'start_write_bigwig',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               MAIN => [
                'bw_seed_all_execution_plans',
                'seed_bigwig_jobs_for_control_experiments'
              ]
            },
        },
        {   -logic_name  => 'bw_seed_all_execution_plans',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::SeedAllExecutionPlans',
            -flow_into   => {
               2 => 'seed_bigwig_jobs_for_signals',
            },
        },
        {   -logic_name  => 'seed_bigwig_jobs_for_signals',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::SeedBigWigJobs',
            -flow_into   => {
               2 => 'write_bigwig_signals',
            },
        },
        {   -logic_name  => 'write_bigwig_signals',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::WriteBigWig',
            -rc_name     => '32Gb_job_3cpus',
            # Allow MEMLIMIT issues to be handled by beekeeper
            -max_retry_count => 1,
            -flow_into   => {
               MAIN     => 'register_signal',
               MEMLIMIT => 'write_bigwig_signal_himem',
            },
        },
        {   -logic_name  => 'write_bigwig_signal_himem',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::WriteBigWig',
            -rc_name     => '64Gb_job_3cpus',
            -flow_into   => {
               MAIN => 'register_signal',
            },
        },
        {   -logic_name  => 'seed_bigwig_jobs_for_control_experiments',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::SeedAllExecutionPlans',
            -flow_into   => {
               2 => 'start_write_bigwig_controls',
               #4 => 'seed_bigwig_jobs_from_list',
            },
        },
        {   -logic_name  => 'start_write_bigwig_controls',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               MAIN => 'seed_bigwig_control_jobs',
            },
        },
        {   -logic_name  => 'seed_bigwig_control_jobs',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::SeedBigWigControlJobs',
            -flow_into   => {
               2 => 'write_control_bigwig',
            },
        },
        {   -logic_name  => 'write_control_bigwig',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::WriteBigWig',
            -rc_name     => '32Gb_job_3cpus',
            # Allow MEMLIMIT issues to be handled by beekeeper
            -max_retry_count => 1,
            -flow_into   => {
               MAIN     => 'register_control_signal',
               MEMLIMIT => 'write_control_bigwig_himem',
            },
        },
        {   -logic_name  => 'write_control_bigwig_himem',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::WriteBigWig',
            -rc_name     => '64Gb_job_3cpus',
            -flow_into   => {
               MAIN => 'register_control_signal',
            },
        },
        {   -logic_name  => 'register_control_signal',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::RegisterSignal',
        },
        {   -logic_name  => 'register_signal',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::RegisterSignal',
        },
    ];
}

1;

