package Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::Base;

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

sub resource_classes {
    my ($self) = @_;
    return {
        %{$self->SUPER::resource_classes},

#         'default' => {
#           'LSF'   => ['-q production-rh7 ', '--reg_conf '.$self->o('reg_conf')], 
#           'LOCAL' => ['', '--reg_conf '.$self->o('reg_conf')] 
#         },
        'default'    => {'LSF' => '-q production-rh7' },
        '250Mb_job'    => {'LSF' => '-q production-rh7 -M250   -R"select[mem>250]   rusage[mem=250]"' },
        '500Mb_job'    => {'LSF' => '-q production-rh7 -M500   -R"select[mem>500]   rusage[mem=500]"' },
        '1Gb_job'      => {'LSF' => '-q production-rh7 -M1000  -R"select[mem>1000]  rusage[mem=1000]"' },
        '2Gb_job'      => {'LSF' => '-q production-rh7 -M2000  -R"select[mem>2000]  rusage[mem=2000]"' },
        '4Gb_job'      => {'LSF' => '-q production-rh7 -M4000  -R"select[mem>4000]  rusage[mem=4000]"' },
        '8Gb_job'      => {'LSF' => '-q production-rh7 -M8000  -R"select[mem>8000]  rusage[mem=8000]"' },
        '16Gb_job'     => {'LSF' => '-q production-rh7 -M16000 -R"select[mem>16000] rusage[mem=16000]"' },
        '24Gb_job'     => {'LSF' => '-q production-rh7 -M24000 -R"select[mem>24000] rusage[mem=24000]"' },
        '32Gb_job'     => {'LSF' => '-q production-rh7 -M32000 -R"select[mem>32000] rusage[mem=32000]"' },
        '48Gb_job'     => {'LSF' => '-q production-rh7 -M48000 -R"select[mem>48000] rusage[mem=48000]"' },
        '64Gb_job'     => {'LSF' => '-q production-rh7 -M64000 -R"select[mem>64000] rusage[mem=64000]"' },
    };
}

1;
