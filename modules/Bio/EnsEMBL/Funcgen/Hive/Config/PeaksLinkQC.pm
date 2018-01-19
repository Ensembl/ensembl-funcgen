=pod 

=head1 NAME

    Bio::EnsEMBL::Funcgen::Hive::Config::PeaksLinkQC

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 LICENSE

    Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
    Copyright [2016-2018] EMBL-European Bioinformatics Institute

    Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

         http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software distributed under the License
    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and limitations under the License.

=head1 CONTACT

=cut
package Bio::EnsEMBL::Funcgen::Hive::Config::PeaksLinkQC;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('Bio::EnsEMBL::Funcgen::Hive::Config::Base');

sub pipeline_analyses {
    my ($self) = @_;
    return [
      
      {
        -logic_name    => 'done_peak_calling',
        -flow_into => {
          MAIN => 'start_qc_proportion_of_reads_in_peaks'
        },
      },
      {   -logic_name => 'start_qc_proportion_of_reads_in_peaks',
	  -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      },
   ];
}

1;
