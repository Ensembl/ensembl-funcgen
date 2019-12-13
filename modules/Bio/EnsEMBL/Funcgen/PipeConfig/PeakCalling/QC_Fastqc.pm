=pod 

=head1 NAME

    Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::QC_Fastqc

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

package Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::QC_Fastqc;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base 'Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::Base';

sub pipeline_analyses {
    my ($self) = @_;
    return [
        {   -logic_name => 'start_fastqc',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into => { 
              'MAIN->A' => [ 'fastqc_seed_signal_experiments', 'fastqc_seed_control_experiments' ],
              'A->MAIN' => [ 'fastqc_done' ],
            },
        },
        {   -logic_name => 'fastqc_seed_signal_experiments',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::SeedAllSignalExperimentNames',
            -flow_into   => {
               2 => 'fastqc_job_factory',
            },
        },
        {   -logic_name => 'fastqc_seed_control_experiments',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::SeedAllControlExperimentNames',
            -flow_into   => {
               2 => 'fastqc_job_factory',
            },
        },
        {   -logic_name => 'fastqc_job_factory',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::QcFastQcJobFactory',
            -parameters => {
              tempdir => '#tempdir_peak_calling#/#species#/fastqc',
            },
            -flow_into => { 
                2 => 'RunFastQC',
            },
        },
        {   -logic_name => 'RunFastQC',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::QcRunFastqc',
            #-rc_name => '2Gb_job',
            # For some reason this analysis is very attractive to workers. 
            # So to prevent them from focusing on this at the beginning of the
            # pipeline, I'm limiting this.
            #
            -analysis_capacity => 5,
            -flow_into => {
                2 => 'QcLoadFastQcResults',
            },
        },
        {   -logic_name        => 'QcLoadFastQcResults',
            -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -analysis_capacity => 1,
            -parameters => {
                use_bash_pipefail => 1,
                cmd => 
                    q( load_fastqc_summary_file.pl             )
                  . q(   --registry      #reg_conf#            )
                  . q(   --species       #species#             )
                  . q(   --read_file_id  #read_file_id#        )
                  . q(   --summary_file  #fastqc_summary_file# )
                  . q(   --work_dir      #fastqc_tempdir#      )
            },
        },
        {   -logic_name => 'fastqc_done',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
    ];
}

1;
