=pod 

=head1 NAME

    Bio::EnsEMBL::Funcgen::Hive::Config::QC_Flagstats

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 LICENSE

    Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');

sub pipeline_analyses {
    my ($self) = @_;
    return [
      {
	-logic_name => 'MergeControlAlignments',
	-flow_into => {
	  MAIN => {
	    BamFileQc => INPUT_PLUS({
		'is_control' => 1,
		'source' => 'MergeControlAlignments',
	      } 
	    ) 
	  },
	},
      },
      {
	-logic_name => 'MergeAlignments',
	-flow_into => {
	  MAIN => { 
	    BamFileQc => INPUT_PLUS({
	      'is_control' => 0,
	      'source' => 'MergeAlignments',
	      } 
	    ) 
	  },
	},
      },
      {
	-logic_name => 'MergeReplicateAlignments',
	-flow_into => {
	  MAIN => { 
	    BamFileQc => INPUT_PLUS({
	      'is_control' => 0,
	      'source' => 'MergeReplicateAlignments',
	      }
	    )
	  },
	},
      },
      {
	-logic_name => 'BamFileQc',
	-module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
	-flow_into => {
	  MAIN => 'QcFlagstatsJobFactory'
	},
	-meadow_type=> 'LOCAL',
      },
      {   -logic_name => 'QcFlagstatsJobFactory',
	  -module     => 'Bio::EnsEMBL::Funcgen::Hive::QcFlagstatsJobFactory',
	  -meadow_type=> 'LOCAL',
	  -flow_into => { 
	    2 => 'QcRunFlagstats',
	  },
      },
      {   -logic_name => 'QcRunFlagstats',
	  -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
	  -meadow_type=> 'LSF',
	  -parameters => { 
		cmd => qq!samtools flagstat #bam_file# > #flagstats_file#!,
	  },
	  -flow_into => { 
	    MAIN => 'LoadFlagstatsToDB',
	  },
      },
      {   -logic_name => 'LoadFlagstatsToDB',
	  -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
	  -meadow_type=> 'LOCAL',
	  -parameters => {
	      cmd => qq(load_samtools_flagstats.pl )
		. qq( --result_set_id #result_set_id# )
		. qq( --flagstats_file #flagstats_file# )
		. qq( --user #tracking_db_user# --pass #tracking_db_pass# --host #tracking_db_host# --dbname #tracking_db_name# )

	  },
      },
    ];
}

1;


