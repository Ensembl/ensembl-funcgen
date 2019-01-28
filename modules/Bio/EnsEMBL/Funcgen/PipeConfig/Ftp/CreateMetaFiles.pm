package Bio::EnsEMBL::Funcgen::PipeConfig::Ftp::CreateMetaFiles;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Funcgen::PipeConfig::Ftp::Base';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub pipeline_analyses {
    my $self = shift;

    return [
        {   -logic_name  => 'start_create_meta_files',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               MAIN => 'create_readme_files'
            },
        },
        {   -logic_name  => 'create_readme_files',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::Ftp::WriteReadmeFile',
            -parameters  => {
                destination => '#ftp_base_dir#/#species#/README',
            },
            -flow_into   => {
               MAIN => 'find_directories_for_checksumming'
            },
        },
        {   -logic_name  => 'find_directories_for_checksumming',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -parameters  => {
                inputcmd     => 'find #ftp_base_dir#/#species# -type d',
                column_names => [ 'directory' ],
            },
            -flow_into   => {
               2 => 'create_checksums'
            },
        },
        {   -logic_name  => 'create_checksums',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => q(cd #directory#; find -maxdepth 1 -type f | grep -v CHECKSUMS | xargs --no-run-if-empty md5sum > CHECKSUMS),
                use_bash_errexit => 1,
                # Can't use this, because of grep in the pipe
                #use_bash_pipefail => 1,
            },
            -flow_into   => {
               MAIN => 'remove_empty_checksum_files'
            },
        },
        {   -logic_name  => 'remove_empty_checksum_files',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                # -f necessary, because there may be no empty checksum files and this
                # prevents rm from giving a non zero exit code.
                #
                cmd => 'find #directory# -size 0b -name CHECKSUMS | xargs rm -f',
                use_bash_errexit => 1,
                use_bash_pipefail => 1,
            },
        },
    ]
}

1;
