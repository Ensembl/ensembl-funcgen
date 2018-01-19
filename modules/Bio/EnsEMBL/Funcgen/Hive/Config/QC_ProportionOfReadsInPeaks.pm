=pod 

=head1 NAME

    Bio::EnsEMBL::Funcgen::Hive::Config::QC_ProportionOfReadsInPeaks

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
package Bio::EnsEMBL::Funcgen::Hive::Config::QC_ProportionOfReadsInPeaks;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('Bio::EnsEMBL::Funcgen::Hive::Config::Base');

sub pipeline_analyses {
    my ($self) = @_;
    return [
      {   -logic_name => 'start_qc_proportion_of_reads_in_peaks',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -flow_into => { 
            'MAIN->A' => 'QcProportionOfReadsInPeaksJobFactory',
            'A->MAIN' => 'done_qc_proportion_of_reads_in_peaks',
          },
      },
      {   -logic_name => 'QcProportionOfReadsInPeaksJobFactory',
	  -module     => 'Bio::EnsEMBL::Funcgen::Hive::QcProportionOfReadsInPeaksJobFactory',
	  -flow_into => { 
	    2 => 'QcProportionOfReadsInPeaks',
	  },
      },
      {   -logic_name => 'QcProportionOfReadsInPeaks',
	  -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
	  -parameters => {
	      cmd => qq( proportion_of_reads_in_peaks.pl )
	      . qq( --peak_file #peak_file#              )
	      . qq( --temp_dir #temp_dir#                )
	      . qq( --peak_caller #peak_caller#          )
	      . qq( --feature_set_id #feature_set_id#    )
	      . qq( --user   #tracking_db_user#   )
              . qq( --pass   #tracking_db_pass#   )
              . qq( --port   #tracking_db_port#   )
	      . qq( --host   #tracking_db_host#   )
	      . qq( --dbname #tracking_db_name#   )
	      . qq( --bam_file #bam_file#         )
	  },
	  -rc_name => '10gb_1cpu',
      },
      {   -logic_name => 'done_qc_proportion_of_reads_in_peaks',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      },
    ];
}

1;
