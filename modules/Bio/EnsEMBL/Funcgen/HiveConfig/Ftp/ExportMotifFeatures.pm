package Bio::EnsEMBL::Funcgen::HiveConfig::Ftp::ExportMotifFeatures;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Hive::PipeConfig::EnsemblGeneric_conf';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
      %{$self->SUPER::pipeline_wide_parameters},

      ftp_base_dir  => $self->o('ftp_base_dir'),
      reg_conf      => $self->o('reg_conf'),
    };
}

sub beekeeper_extra_cmdline_options {
    my ($self) = @_;
    return '-reg_conf ' . $self->o('reg_conf') . ' -keep_alive -can_respecialize 1';
}

sub pipeline_analyses {
    my $self = shift;

    return [
        {   -logic_name  => 'start_export',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               'MAIN' => 'export_motif_features'
            },
        },
        {   -logic_name  => 'export_motif_features',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'export_motif_features.pl --ftp_base_dir #ftp_base_dir#/#species# --registry #reg_conf# --species #species#',

            },
        },
    ]
}

1;
