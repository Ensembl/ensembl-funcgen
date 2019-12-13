=pod 

=head1 NAME

    Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::ChIPSeqCleanup

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 LICENSE

    Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
    Copyright [2016-2020] EMBL-European Bioinformatics Institute

    Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

         http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software distributed under the License
    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and limitations under the License.

=head1 CONTACT

=cut

package Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::ChIPSeqCleanup;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base 'Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::Base';

sub pipeline_analyses {
    my ($self) = @_;
    return [
        {   -logic_name => 'start_cleanup',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into => { 
              MAIN => 'remove_intermediary_bam_files',
            },
        },
        {   -logic_name => 'remove_intermediary_bam_files',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd => qq(remove_intermediary_bam_files.pl --registry #reg_conf# --species #species# --data_root_dir #data_root_dir#)
            },
            -flow_into => {
                MAIN => 'remove_unregistered_files',
            },
        },
        {   -logic_name => 'remove_unregistered_files',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd => qq(remove_unregistered_files.pl --registry #reg_conf# --species #species# --data_root_dir #data_root_dir#)
            },
            -flow_into => {
                MAIN => 'remove_empty_directories',
            },
        },
        {   -logic_name => 'remove_empty_directories',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd => qq(remove_empty_directories.pl --registry #reg_conf# --species #species# --data_root_dir #data_root_dir#)
            },
            -flow_into => {
                MAIN => 'cleanup_done',
            },
        },
        {   -logic_name => 'cleanup_done',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
    ];
}

1;
