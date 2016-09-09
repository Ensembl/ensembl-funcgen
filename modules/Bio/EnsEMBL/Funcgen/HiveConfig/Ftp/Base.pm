package Bio::EnsEMBL::Funcgen::HiveConfig::Ftp::Base;

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
      temp_dir      => '#ftp_base_dir#/#species#/tempdir',
    };
}

sub beekeeper_extra_cmdline_options {
    my ($self) = @_;
    return '-reg_conf ' . $self->o('reg_conf') . ' -keep_alive -can_respecialize 1';
}

1;
