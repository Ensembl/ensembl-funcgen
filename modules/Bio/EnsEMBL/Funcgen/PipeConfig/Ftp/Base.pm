package Bio::EnsEMBL::Funcgen::PipeConfig::Ftp::Base;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::ResourceClasses';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub default_options {
    my ($self) = @_;
    
    use Date::Format;
    # Looks like: 20160928
    my $data_freeze_date = time2str('%Y%m%d', time);

    return {
        %{$self->SUPER::default_options},
        data_freeze_date => $data_freeze_date,
        tempdir => '',
    }
}

sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
      %{$self->SUPER::pipeline_wide_parameters},

      ftp_base_dir     => $self->o('ftp_base_dir'),
      reg_conf         => $self->o('reg_conf'),
      data_freeze_date => $self->o('data_freeze_date'),
      tempdir          => $self->o('tempdir'),
      tempdir_ftp      => $self->o('tempdir') . '/ftp_export',
    };
}

sub beekeeper_extra_cmdline_options {
    my ($self) = @_;
    return '-reg_conf ' . $self->o('reg_conf') . ' -keep_alive -can_respecialize 1';
}

1;
