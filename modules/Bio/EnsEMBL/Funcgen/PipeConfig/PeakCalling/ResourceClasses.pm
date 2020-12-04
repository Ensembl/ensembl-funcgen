package Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::ResourceClasses;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub _gigabytes_to_bytes {

    my $self      = shift;
    my $gigabytes = shift;
    
    return $gigabytes * 1073741824;
}

sub _cpu_to_nanocpu {

    my $self = shift;
    my $cpu = shift;
    
    return $cpu * 1000000000;
}

sub _swarm_resource {

    my $self = shift;

    my $cpus = shift;
    my $gb   = shift;
    
    use Bio::EnsEMBL::Hive::Utils 'stringify';
    
    my $gb_value  = $self->_gigabytes_to_bytes ( $gb   );
    my $cpu_value = $self->_cpu_to_nanocpu     ( $cpus );

    return stringify(
        {
            'Limits'  => {
                'NanoCPUs'     => $self->_cpu_to_nanocpu     ( $cpus ),
                'MemoryBytes'  => $self->_gigabytes_to_bytes ( $gb   ),
            },
            'Reservations'  => {
                'NanoCPUs'     => $self->_cpu_to_nanocpu     ( $cpus ),
                'MemoryBytes'  => $self->_gigabytes_to_bytes ( $gb   ),
            },
        }
    )
}

sub resource_classes {
    my ($self) = @_;

    my $lsf_queue_option = '-q production-rh74';
    
    my $resource_classes = {
        %{$self->SUPER::resource_classes},  # inherit 'default' from the parent class

        'default'         => {'LSF' => $lsf_queue_option   },
        '250Mb_job'       => {
            'LSF'         => [ qq($lsf_queue_option -M250   -R"select[mem>250]   rusage[mem=250]") ],
            'DockerSwarm' => [ $self->_swarm_resource(1, 0.25) ],
        },
        '500Mb_job'       => {
            'LSF'         => [ qq($lsf_queue_option -M500   -R"select[mem>500]   rusage[mem=500]") ],
            'DockerSwarm' => [ $self->_swarm_resource(1, 0.5) ],
        },
        '1Gb_job'         => {
            'LSF'         => [ qq($lsf_queue_option -M1000  -R"select[mem>1000]  rusage[mem=1000]") ],
            'DockerSwarm' => [ $self->_swarm_resource(1, 1 - 0.5) ],
        },
        '2Gb_job'         => {
            'LSF'         => [ qq($lsf_queue_option -M2000  -R"select[mem>2000]  rusage[mem=2000]") ],
            'DockerSwarm' => [ $self->_swarm_resource(1, 2 -1) ],
        },
        '4Gb_job'         => {
            'LSF'         => [ qq($lsf_queue_option -M4000  -R"select[mem>4000]  rusage[mem=4000]") ],
            'DockerSwarm' => [ $self->_swarm_resource(1, 4 -1) ],
        },
        '4Gb_job_2h'      => {
            'LSF'         => [ qq($lsf_queue_option -M4000  -R"select[mem>4000]  rusage[mem=4000]" -W 2:00) ],
            'DockerSwarm' => [ $self->_swarm_resource(1, 4 -1) ],
        },
        '4Gb_job_2cpus'   => {
            'LSF'         => qq($lsf_queue_option -M4000  -R"select[mem>4000]  rusage[mem=4000]" -n 2),
            'DockerSwarm' => [ $self->_swarm_resource(1, 4 -1) ],
        },
        '8Gb_job'         => {
            'LSF'         => [ qq($lsf_queue_option -M8000  -R"select[mem>8000]  rusage[mem=8000]") ],
            'DockerSwarm' => [ $self->_swarm_resource(1, 8 -1) ],
        },
        '8Gb_job_2h'     => {
            'LSF'         => qq($lsf_queue_option -M8000 -R"select[mem>8000] rusage[mem=8000]" -W 2:00),
            'DockerSwarm' => [ $self->_swarm_resource(1, 8 -1) ],
        },
        '8Gb_job_8h'     => {
            'LSF'         => qq($lsf_queue_option -M8000 -R"select[mem>8000] rusage[mem=8000]" -W 8:00),
            'DockerSwarm' => [ $self->_swarm_resource(1, 8 -1) ],
        },
        '10Gb_job_2cpus'   => {
            'LSF'         => qq($lsf_queue_option -M10000  -R"select[mem>10000]  rusage[mem=10000]" -n 2),
            'DockerSwarm' => [ $self->_swarm_resource(1, 10 -1) ],
        },
	'16Gb_job'        => {
            'LSF'         => [ qq($lsf_queue_option -M16000 -R"select[mem>16000] rusage[mem=16000]") ],
            'DockerSwarm' => [ $self->_swarm_resource(1, 16 -1) ],
        },
        '24Gb_job'        => {
            'LSF'         => [ qq($lsf_queue_option -M24000 -R"select[mem>24000] rusage[mem=24000]") ],
            'DockerSwarm' => [ $self->_swarm_resource(1, 24 -1) ],
        },
        '32Gb_job'        => {
            'LSF'         => [ qq($lsf_queue_option -M32000 -R"select[mem>32000] rusage[mem=32000]") ],
            'DockerSwarm' => [ $self->_swarm_resource(1, 32 -1) ],
        },
        '32Gb_job_2h'     => {
            'LSF'         => qq($lsf_queue_option -M32000 -R"select[mem>32000] rusage[mem=32000]" -W 2:00),
            'DockerSwarm' => [ $self->_swarm_resource(1, 32 -1) ],
        },
        '32Gb_job_8h'     => {
            'LSF'         => qq($lsf_queue_option -M32000 -R"select[mem>32000] rusage[mem=32000]" -W 8:00),
            'DockerSwarm' => [ $self->_swarm_resource(1, 32 -1) ],
        },
        '32Gb_job_2cpus'  => {
            'LSF'         => qq($lsf_queue_option -M32000 -R"select[mem>32000] rusage[mem=32000]" -n 2),
            'DockerSwarm' => [ $self->_swarm_resource(2, 32 -1) ],
        },
        '32Gb_job_3cpus'  => {
            'LSF'         => qq($lsf_queue_option -M32000 -R"select[mem>32000] rusage[mem=32000]" -n 3),
            'DockerSwarm' => [ $self->_swarm_resource(3, 32 -1) ],
        },
        '48Gb_job'        => {
            'LSF'         => [ qq($lsf_queue_option -M48000 -R"select[mem>48000] rusage[mem=48000]") ],
            'DockerSwarm' => [ $self->_swarm_resource(1, 48 -1) ],
        },
        '64Gb_job'        => {
            'LSF'         => [ qq($lsf_queue_option -M64000 -R"select[mem>64000] rusage[mem=64000]") ],
            'DockerSwarm' => [ $self->_swarm_resource(1, 64 -1 ) ],
        },
        '64Gb_job_3cpus'  => {
            'LSF'         => qq($lsf_queue_option -M64000 -R"select[mem>64000] rusage[mem=64000]" -n 3),
            'DockerSwarm' => [ $self->_swarm_resource(3, 64 -1 ) ],
        },
        'binarization'    => {
            'LSF'         => qq($lsf_queue_option -M16000 -R"select[mem>16000] rusage[mem=16000]" -n 5),
            'DockerSwarm' => [ $self->_swarm_resource(5, 16 -1 ) ],
        },
        'learn_model'     => {
            'LSF'         => qq($lsf_queue_option -M31000 -R"span[hosts=1] select[mem>31000] rusage[mem=31000]" -n 12),
            'DockerSwarm' => [ $self->_swarm_resource(12, 31 -1 ) ],
        },
    };
    return $resource_classes;
}

1;
