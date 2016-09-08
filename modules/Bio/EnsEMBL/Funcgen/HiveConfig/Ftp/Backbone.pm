package Bio::EnsEMBL::Funcgen::HiveConfig::Ftp::Backbone;

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
      temp_dir      => '#ftp_base_dir#/tempdir/annotated_features',
    };
}

sub beekeeper_extra_cmdline_options {
    my ($self) = @_;
    return '-reg_conf ' . $self->o('reg_conf') . ' -keep_alive -can_respecialize 1';
}

sub pipeline_analyses {
    my $self = shift;

    return [
        {   -logic_name  => 'backbone_fire_exports',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_export',
               'A->1' => 'backbone_fire_createmetafiles'
            },
        },
        {   -logic_name  => 'backbone_fire_createmetafiles',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'create_readme_files',
               'A->1' => 'backbone_pipeline_finished'
            },
        },
        {   -logic_name => 'backbone_pipeline_finished',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        }
    ]
}

1;
