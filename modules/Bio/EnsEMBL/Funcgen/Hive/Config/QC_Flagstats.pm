=pod 

=head1 NAME

    Bio::EnsEMBL::Funcgen::Hive::Config::QC_Flagstats

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
package Bio::EnsEMBL::Funcgen::Hive::Config::QC_Flagstats;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('Bio::EnsEMBL::Funcgen::Hive::Config::Base');

sub pipeline_analyses {
    my ($self) = @_;
    return [
      {
	-logic_name => 'BamFileQc',
	-flow_into => {
	  MAIN => WHEN(
	    '#has_unmapped_reads# eq "yes"' => 'qc_flagstats_start',
	  ),
	},
      },
        {   -logic_name => 'qc_flagstats_start',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into => { 
              'MAIN->A' => 'QcFlagstatsJobFactory',
              'A->MAIN' => 'qc_flagstats_done',
            },
        },
      {   -logic_name => 'QcFlagstatsJobFactory',
	  -module     => 'Bio::EnsEMBL::Funcgen::Hive::QcFlagstatsJobFactory',
	  -flow_into => { 
	    2 => 'QcRunFlagstats',
	  },
      },
      {   -logic_name => 'QcRunFlagstats',
	  -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
	  -parameters => { 
		cmd => qq!samtools flagstat #bam_file# > #flagstats_file#!,
	  },
	  -flow_into => { 
	    MAIN => 'LoadFlagstatsToDB',
	  },
      },
      {   -logic_name => 'LoadFlagstatsToDB',
	  -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
	  -parameters => {
	      cmd => qq(load_samtools_flagstats.pl )
		. qq( --result_set_id #result_set_id# )
		. qq( --flagstats_file #flagstats_file# )
		. qq( --user #tracking_db_user# --pass #tracking_db_pass# --host #tracking_db_host# --port #tracking_db_port# --dbname #tracking_db_name# )
		. qq( --work_dir #tempdir#  )
		. qq( --bam_file #bam_file# )
	  },
      },
      {   -logic_name => 'qc_flagstats_done',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      },

#       {   -logic_name => 'DeleteBamWithDuplicates',
# 	  -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
# 	  -parameters => { 
# 		cmd => qq(rm #bam_file#),
# 	  },
#       },

    ];
}

1;
