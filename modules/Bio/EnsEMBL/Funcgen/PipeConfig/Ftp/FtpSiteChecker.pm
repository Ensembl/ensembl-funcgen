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
              MAIN => [
                'ftp_check_regulatory_activity_directories_exist',
                'ftp_check_regulatory_activity_files_exist',
                'ftp_check_regulatory_feature_numbers',
                'ftp_check_file_content_consistent_job_factory',
                'ftp_check_activity_summaries_consistent',
              ]
            },
        },
        {   -logic_name  => 'ftp_check_regulatory_activity_directories_exist',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -rc_name     => '32Gb_job',
            -parameters  => {
                cmd => q(
                  ftp_site_checks.pl \
                    --registry #reg_conf# \
                    --species #species# \
                    --ftp_dir #ftp_base_dir#/#species# \
                    --check check_regulatory_activity_directories_exist
                ),
            },
        },
        {   -logic_name  => 'ftp_check_regulatory_activity_files_exist',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -rc_name     => '32Gb_job',
            -parameters  => {
                cmd => q(
                  ftp_site_checks.pl \
                    --registry #reg_conf# \
                    --species #species# \
                    --ftp_dir #ftp_base_dir#/#species# \
                    --check check_regulatory_activity_files_exist
                ),
            },
        },
        {   -logic_name  => 'ftp_check_regulatory_feature_numbers',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -rc_name     => '32Gb_job',
            -parameters  => {
                cmd => q(
                  ftp_site_checks.pl \
                    --registry #reg_conf# \
                    --species #species# \
                    --ftp_dir #ftp_base_dir#/#species# \
                    --check check_regulatory_feature_numbers
                ),
            },
        },
        {   -logic_name  => 'ftp_check_file_content_consistent_job_factory',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::Ftp::JobFactoryFtpRegulatoryActivityFiles',
            -rc_name     => '32Gb_job',
            -flow_into => { 
                2 => 'ftp_check_file_content_consistent',
              },
        },
        {   -logic_name  => 'ftp_check_file_content_consistent',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -rc_name     => '32Gb_job',
            -analysis_capacity => 50,
            -parameters  => {
                cmd => q(
                  ftp_site_checks.pl \
                    --registry #reg_conf# \
                    --species #species# \
                    --ftp_dir #ftp_base_dir#/#species# \
                    --check check_file_content_consistent \
                    --check_file #file_name#
                ),
            },
        },
        {   -logic_name  => 'ftp_check_activity_summaries_consistent',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -rc_name     => '32Gb_job',
            -parameters  => {
                cmd => q(
                  ftp_site_checks.pl \
                    --registry #reg_conf# \
                    --species #species# \
                    --ftp_dir #ftp_base_dir#/#species# \
                    --check check_activity_summaries_consistent
                ),
            },
        },
    ]
}

1;
