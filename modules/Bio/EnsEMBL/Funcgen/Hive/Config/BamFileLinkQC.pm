=pod 

=head1 NAME

    Bio::EnsEMBL::Funcgen::Hive::Config::BamFileLinkQC

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
package Bio::EnsEMBL::Funcgen::Hive::Config::BamFileLinkQC;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('Bio::EnsEMBL::Funcgen::Hive::Config::Base');

sub pipeline_analyses {
    my ($self) = @_;
    return [
      {
	-logic_name => 'done_align_controls',
	-flow_into => {
	  MAIN => {
	    BamFileQc => INPUT_PLUS({
		'is_control' => 1,
		'source' => 'done_align_controls',
		'has_duplicates' => 'yes',
		'has_unmapped_reads' => 'yes',
		'bam_file_for_qc' => '#bam_file_with_unmapped_reads_and_duplicates#',
	      } 
	    ) 
	  },
	},
      },
      {
	-logic_name => 'RemoveDuplicateControlAlignments',
	-flow_into => {
	  MAIN => {
	    BamFileQc => INPUT_PLUS({
		'is_control' => 1,
		'source' => 'RemoveDuplicateControlAlignments',
		'has_duplicates' => 'no',
		'has_unmapped_reads' => 'no',
		'bam_file_for_qc' => '#bam_file#',
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
	      'has_duplicates' => 'yes',
	      'has_unmapped_reads' => 'yes',
	      'bam_file_for_qc' => '#bam_file_with_unmapped_reads_and_duplicates#',
	      } 
	    ) 
	  },
	},
      },
      {
	-logic_name => 'RemoveDuplicateAlignments',
	-flow_into => {
	  MAIN => { 
	    BamFileQc => INPUT_PLUS({
	      'is_control' => 0,
	      'source' => 'RemoveDuplicateAlignments',
	      'has_duplicates' => 'no',
	      'has_unmapped_reads' => 'no',
	      'bam_file_for_qc' => '#bam_file#',
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
	      'has_duplicates' => 'yes',
	      'has_unmapped_reads' => 'yes',
	      'bam_file_for_qc' => '#bam_file_with_unmapped_reads_and_duplicates#',
	      }
	    )
	  },
	},
      },
      {
	-logic_name => 'RemoveDuplicateReplicateAlignments',
	-flow_into => {
	  MAIN => { 
	    BamFileQc => INPUT_PLUS({
	      'is_control' => 0,
	      'source' => 'RemoveDuplicateReplicateAlignments',
	      'has_duplicates' => 'no',
	      'has_unmapped_reads' => 'no',
	      'bam_file_for_qc' => '#bam_file#',
	      }
	    )
	  },
	},
      },
      {
	-logic_name => 'FixReplicateResultSetsExperimentIds',
	-flow_into => {
	  MAIN => { 
	    BamFileQc => INPUT_PLUS({
	      'is_control' => 0,
	      'source' => 'FixReplicateResultSetsExperimentIds',
	      'has_duplicates' => undef,
	      'has_unmapped_reads' => 'no',
	      'bam_file_for_qc' => '#bam_file#',
	      }
	    )
	  },
	},
      },
      {
	-logic_name => 'BamFileQc',
	-module     => 'Bio::EnsEMBL::Funcgen::Hive::BamFileQc',
# 	-flow_into => {
# 	    2 => WHEN(
#                 '#has_duplicates# eq "yes"' => { 
#                   ':////accu?file_to_delete=[]' => { 
#                     'file_to_delete' => '#bam_file#'
#                   } 
#                 }
# 	      ),
# 	  },

      },
    ];
}

1;



