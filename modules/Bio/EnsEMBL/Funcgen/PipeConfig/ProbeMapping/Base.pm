package Bio::EnsEMBL::Funcgen::PipeConfig::ProbeMapping::Base;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

sub pipeline_wide_parameters {
    my $self = shift;
    return {
      %{$self->SUPER::pipeline_wide_parameters},

      tempdir                 => $self->o('tempdir'),
      probe_directory         => $self->o('probe_directory'),
      reg_conf                => $self->o('reg_conf'),
      unmapped_sequences_file => '#tempdir#/#species#/unmapped_probe_sequences.fasta',
      toplevel_sequences_file => '#tempdir#/#species#/toplevel.fasta',
      gene_sequences_file     => '#tempdir#/#species#/genes.fasta',
    };
}

sub default_options {
    my ($self) = @_;
    return {
      %{ $self->SUPER::default_options() },

      tempdir          => '/nfs/nobackup/ensembl/mnuhn/array_mapping/temp',
      probe_directory  => '/nfs/production/panda/ensembl/funcgen/array_mapping/',
      pipeline_name    => 'probemapping',
    };
}

sub beekeeper_extra_cmdline_options {
    my ($self) = @_;
    return '-reg_conf ' . $self->o('reg_conf') . ' -keep_alive -can_respecialize 1';
}

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

    #return qq({"Limits":{"MemoryBytes":$gb_value,"NanoCPUs":$cpu_value},"Reservations":{"NanoCPUs":$cpu_value,"MemoryBytes":$gb_value}});

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
    
    return {
        %{$self->SUPER::resource_classes},

        'default'         => {'LSF' => $lsf_queue_option },
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
        '8Gb_job'         => {
            'LSF'         => [ qq($lsf_queue_option -M8000  -R"select[mem>8000]  rusage[mem=8000]") ],
            'DockerSwarm' => [ $self->_swarm_resource(1, 8 -1) ],
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
        '48Gb_job'        => {
            'LSF'         => [ qq($lsf_queue_option -M48000 -R"select[mem>48000] rusage[mem=48000]") ],
            'DockerSwarm' => [ $self->_swarm_resource(1, 48 -1) ],
        },
        '64Gb_job'        => {
            'LSF'         => [ qq($lsf_queue_option -M64000 -R"select[mem>64000] rusage[mem=64000]") ],
            'DockerSwarm' => [ $self->_swarm_resource(1, 64 -1 ) ],
        },
        'parallel_sort'   => { 
            'LSF'         => [ qq($lsf_queue_option -n 16 -M32000 -R"select[mem>32000] rusage[mem=32000]") ],
            'DockerSwarm' => [ $self->_swarm_resource(16, 32) ],
        },
    };
}

1;
