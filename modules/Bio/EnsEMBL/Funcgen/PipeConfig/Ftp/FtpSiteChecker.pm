package Bio::EnsEMBL::Funcgen::PipeConfig::Ftp::FtpSiteChecker;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Funcgen::PipeConfig::Ftp::Base';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub pipeline_analyses {
    my $self = shift;

    return [
        {   -logic_name => 'start_ftp_site_checker',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into => { 
              MAIN => 'ftp_site_checks',
            },
        },
        {   -logic_name  => 'ftp_site_checks',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => qq(
                  ftp_site_checks.pl \
                    --registry #registry# \
                    --species #species# \
                    --ftp_dir #ftp_base_dir#/#species#
                ),
            },
        },
    ]
}

1;
