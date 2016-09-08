package Bio::EnsEMBL::Funcgen::HiveConfig::Ftp::CreateMetaFiles;

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
        {   -logic_name  => 'create_readme_files',
            -module      => 'Bio::EnsEMBL::Funcgen::Hive::Ftp::WriteReadmeFile',
            -parameters  => {
                destination => '#ftp_base_dir#/#species#/README',
            },
            -flow_into   => {
               'MAIN' => 'create_checksums'
            },
        },
        {   -logic_name  => 'create_checksums',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'find #ftp_base_dir#/#species# -type d | xargs -L 1 -Ipp echo "cd pp; md5sum * > CHECKSUMS" | bash',
            },
            -flow_into   => {
               'MAIN' => 'remove_empty_checksum_files'
            },
        },
        {   -logic_name  => 'remove_empty_checksum_files',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
		# -f necessary, because there may be no empty checksum files and this
		# prevents rm from giving a non zero exit code.
		#
                cmd => 'find #ftp_base_dir#/#species# -size 0b -name CHECKSUMS | xargs rm -f',
            },
        },
    ]
}

1;
