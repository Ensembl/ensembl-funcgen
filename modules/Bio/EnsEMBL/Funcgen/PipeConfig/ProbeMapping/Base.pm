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
    
    my $slurm_queue_option = '--partition=production';
    
    return {
        %{$self->SUPER::resource_classes},

        'default'         => {'SLURM' => qq($slurm_queue_option --mem=250 --time=72:00:00) },
        '250Mb_job'       => {
            'SLURM'         => [ qq($slurm_queue_option --mem=250i --time=72:00:00) ],
            'DockerSwarm' => [ $self->_swarm_resource(1, 0.25) ],
        },
        '500Mb_job'       => {
            'SLURM'         => [ qq($slurm_queue_option --mem=500i --time=72:00:00) ],
            'DockerSwarm' => [ $self->_swarm_resource(1, 0.5) ],
        },
        '1Gb_job'         => {
            'SLURM'         => [ qq($slurm_queue_option --mem=1Gi --time=72:00:00) ],
            'DockerSwarm' => [ $self->_swarm_resource(1, 1 - 0.5) ],
        },
        '2Gb_job'         => {
            'SLURM'         => [ qq($slurm_queue_option --mem=2G --time=72:00:00) ],
            'DockerSwarm' => [ $self->_swarm_resource(1, 2 -1) ],
        },
        '4Gb_job'         => {
            'SLURM'         => [ qq($slurm_queue_option --mem=4G --time=72:00:00) ],
            'DockerSwarm' => [ $self->_swarm_resource(1, 4 -1) ],
        },
        '8Gb_job'         => {
            'SLURM'         => [ qq($slurm_queue_option --mem=8G --time=72:00:00) ],
            'DockerSwarm' => [ $self->_swarm_resource(1, 8 -1) ],
        },
        '16Gb_job'        => {
            'SLURM'         => [ qq($slurm_queue_option --mem=16G --time=72:00:00) ],
            'DockerSwarm' => [ $self->_swarm_resource(1, 16 -1) ],
        },
        '24Gb_job'        => {
            'SLURM'         => [ qq($slurm_queue_option --mem=24G --time=72:00:00) ],
            'DockerSwarm' => [ $self->_swarm_resource(1, 24 -1) ],
        },
        '32Gb_job'        => {
            'SLURM'         => [ qq($slurm_queue_option --mem=32G --time=72:00:00) ],
            'DockerSwarm' => [ $self->_swarm_resource(1, 32 -1) ],
        },
        '48Gb_job'        => {
            'SLURM'         => [ qq($slurm_queue_option --mem=48G --time=72:00:00) ],
            'DockerSwarm' => [ $self->_swarm_resource(1, 48 -1) ],
        },
        '64Gb_job'        => {
            'SLURM'         => [ qq($slurm_queue_option --mem=64G --time=72:00:00) ],
            'DockerSwarm' => [ $self->_swarm_resource(1, 64 -1 ) ],
        },
        'parallel_sort'   => { 
            'SLURM'         => [ qq($slurm_queue_option -n 16 --mem=32G --time=72:00:00) ],
            'DockerSwarm' => [ $self->_swarm_resource(16, 32) ],
        },
    };
}

1;
