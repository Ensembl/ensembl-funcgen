=pod 

=head1 NAME

    Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::QC_Reports

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 LICENSE

    Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
    Copyright [2016-2022] EMBL-European Bioinformatics Institute

    Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

         http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software distributed under the License
    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and limitations under the License.

=head1 CONTACT

=cut
package Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::QC_FastqcReport;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base 'Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::Base';

sub pipeline_analyses {
    my ($self) = @_;
    return [
        {   -logic_name => 'start_fastqc_report',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into => { 
              '1->A' => 'generate_fastqc_report',
              'A->1' => 'fastqc_report_done',
            },
        },
        {   -logic_name  => 'generate_fastqc_report',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                  cmd => 
                      q(
                        generate_fastqc_report.pl \
                            --registry #reg_conf# \
                            --species #species# \
                            --output_file #reports_dir#/#species#/fastqc.html
                      )
            },
          -rc_name    => '1Gb_job',
        },
        {   -logic_name => 'fastqc_report_done',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
    ];
}

1;
