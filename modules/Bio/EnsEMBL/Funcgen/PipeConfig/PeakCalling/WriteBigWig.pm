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
               MAIN => 'seed_bigwig_jobs',
            },
        },
        {   -logic_name  => 'seed_bigwig_jobs',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::SeedBigWigJobs',
            -flow_into   => {
               2 => 'write_bigwig',
            },
        },
        {   -logic_name  => 'write_bigwig',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::WriteBigWig',
            -rc_name     => '32Gb_job_3cpus',
            # Allow MEMLIMIT issues to be handled by beekeeper
            -max_retry_count => 1,
            -flow_into   => {
               MAIN     => 'register_signal',
               MEMLIMIT => 'write_bigwig_himem',
            },
        },
        {   -logic_name  => 'write_bigwig_himem',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::WriteBigWig',
            -rc_name     => '64Gb_job_3cpus',
            -flow_into   => {
               MAIN => 'register_signal',
            },
        },
        {   -logic_name  => 'register_signal',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::RegisterSignal',
        },
    ];
}

1;

