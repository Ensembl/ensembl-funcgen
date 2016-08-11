package Bio::EnsEMBL::Funcgen::Hive::Config::ReadAlignment;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('Bio::EnsEMBL::Funcgen::Hive::Config::BaseSequenceAnalysis');

sub default_options {
  my $self = shift;
  return {
    %{$self->SUPER::default_options},

      #Size of each sequence chunk to be aligned (nbr of reads * 4)
      #fastq_chunk_size      => 16000000, #This should run in 30min-1h
      fastq_chunk_size      =>     5000000,
#       fastq_chunk_size      =>  20000000,
      alignment_analysis    => 'bwa_samse',
   };
}

sub pipeline_wide_parameters {
    my $self = shift;
    return {
      %{$self->SUPER::pipeline_wide_parameters},

      #Size of each sequence chunk to be aligned (nbr of reads * 4)
      fastq_chunk_size      => $self->o('fastq_chunk_size'),
      alignment_analysis    => $self->o('alignment_analysis'),
    };
}

sub pipeline_analyses {
  my $self = shift;
  return
   [
    @{$self->SUPER::pipeline_analyses}, #To pick up BaseSequenceAnalysis-DefineMergedOutputSet

    {
	-logic_name => 'PrePipelineChecks',
	-input_ids  => [ {} ],
	-module     => 'Bio::EnsEMBL::Funcgen::Hive::ErsaPrePipelineChecks',
    },
    {   -logic_name => 'JobPool',
	-module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
	-wait_for    => 'JobPool',
	-flow_into => {
	  MAIN => 'TokenLimitedJobFactory',
	},
    },
    {   -logic_name => 'TokenLimitedJobFactory',
	-module     => 'Bio::EnsEMBL::Funcgen::Hive::TokenLimitedJobFactory',
	-flow_into => {
	  '2->A' => 'Dummy',
	  'A->1' => 'TokenLimitedJobFactory',
	},
    },
    {
      -logic_name => 'Dummy',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
	-flow_into => {
	  MAIN => WHEN(
	    '#identify_controls# == 1'  => { 'IdentifyAlignInputSubsets' => INPUT_PLUS() },
	    ELSE { 'CreateJobBatchUsingNewGroupingMechanism' => INPUT_PLUS() },
	  ),
	},
    },
    {
      -logic_name => 'CreateJobBatchUsingNewGroupingMechanism',
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::CreateJobBatchUsingNewGroupingMechanism',
      -parameters => {
      },
      -wait_for => 'PrePipelineChecks',
      -flow_into => {
	'2->A' => 'DefineResultSets',
	'A->2' => 'DeleteFilesFromJobFan',
      },
    },
    {
      -logic_name => 'IdentifyAlignInputSubsets',
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::IdentifySetInputs',
      -parameters => {
      set_type                     => 'InputSubset',
	feature_set_analysis_type    => 'peak',
	default_feature_set_analyses => $self->o('default_peak_analyses'),
	dataflow_param_names => ['no_idr'], 
      },
      -wait_for => 'PrePipelineChecks',
      -flow_into => {
	'2->A' => 'DefineResultSets',
	'A->4' => 'DeleteFilesFromJobFan',
      },
    },
    {
      -logic_name => 'DeleteFilesFromJobFan',
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::ErsaCleanup',
    },
    {
     -logic_name => 'DefineResultSets',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::DefineResultSets',
    -flow_into => {
      '2->A' => 'FixResultSetsExperimentIds',
      'A->3' => 'SplitFastqFilesFromControls',
     },
    },
    {
     -logic_name => 'FixResultSetsExperimentIds',
     -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
     -parameters => {
	#
	# Sets the experiment id for the current result set.
	#
	# This should be set when creating the result set, but it is not. 
	# Until this is fixed in the api we do it here in an extra step.
	#
	# HACK the "and result_set.name like concat(experiment.name, "%")" bit is necessary, because there can be more than on experiment with the same feature type and epigenome, so to meet deadline using the fact that the names match up by convention. But this has to be set properly.
	#
	sql => qq(
	  update result_set, epigenome, experiment 
	  set result_set.experiment_id = experiment.experiment_id
	  where epigenome.epigenome_id=experiment.epigenome_id and epigenome.epigenome_id=result_set.epigenome_id and experiment.feature_type_id=result_set.feature_type_id and result_set.name like concat(experiment.name, "%")
	  and result_set_id = #dbID#
	  ),
	  db_conn => '#out_db_url#'
	},
#       -meadow     => 'LOCAL',
      -analysis_capacity => 1,
    },
    {
      -logic_name => 'SplitFastqFilesFromControls',
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::PreprocessFastqs',
      -flow_into => {
	'2->A' => 'AlignChunksFromControls',
	'A->3' => 'MergeControlAlignments',
	},
      # BWA only uses one cpu, but we are asking for two.
      #
      # The reason is that in extremely rare cases jobs were suspended by 
      # LSF. This happened in 3 out of 20000 jobs, in one test run. Often this
      # doesn't happen at all. However, when this happens, the pipeline can't
      # move on, so we are asking for an extra cpu for these extremely rare
      # cases.
      #
      -rc_name => '10gb_2cpu'
     },
     {
      -logic_name => 'SplitMergedFastQ',
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::PreprocessFastqs',
      -parameters => { 
	merge => 1 
      },
      -flow_into => {
	'2->A' => 'AlignChunksFromMergedFastqs',
	'A->3' =>  'MergeAlignments',
	},
      -rc_name => '10gb_2cpu'
     },
     {
      -logic_name => 'SplitFastqsFromReplicatedExperiments',
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::PreprocessFastqs',
      -flow_into => {
	'2->A' => 'AlignChunksFromReplicateExperiments',
	'A->3' => 'MergeReplicateAlignments' ,
	},
      -rc_name => '10gb_2cpu'
     },
    {
      -logic_name => 'AlignChunksFromControls',
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::RunAligner',
      -rc_name => 'normal_10gb'
     },
    {
    -logic_name => 'AlignChunksFromMergedFastqs',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::RunAligner',
     -rc_name => 'normal_10gb'
     },
    {
      -logic_name => 'AlignChunksFromReplicateExperiments',
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::RunAligner',
      -rc_name => 'normal_10gb'
     },
    {
      -logic_name => 'MergeControlAlignments',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::MergeAlignments',
     -parameters => {
	run_controls => 1,
     },
     -flow_into => {
	  MAIN => {
	    'RemoveDuplicateControlAlignments' => undef,
	    ':////accu?file_to_delete=[]' => { 
	      'file_to_delete' => '#bam_file_with_unmapped_reads_and_duplicates#'
	    },
# 	    # Create bigwigs for controls
#  	    'WriteBigWig' => undef,
	  },
       },
     -rc_name => 'normal_monitored_2GB',
    },
    {
     -logic_name => 'RemoveDuplicateControlAlignments',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::RemoveDuplicateAlignments',
     -parameters => {
	run_controls => 1,
     },
     -flow_into => {
	 MAIN => [ 
	  'JobFactorySignalProcessing', 
	  # Create bigwigs for controls
	  WHEN('1' => { 'WriteBigWig' => INPUT_PLUS({ type => 'control' })})
         ],
       },
#      -rc_name => '64GB_3cpu',
      -rc_name => 'normal_4GB_2cpu',
    },
    {
      -logic_name => 'JobFactorySignalProcessing',
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::JobFactorySignalProcessing',
      -flow_into => {
	'A->3'  => 'CleanupFilesFromPermissiveSWEmblJobFan',
	'10'    => 'SplitMergedFastQ' ,
	'11->A' => 'SplitFastqsFromReplicatedExperiments',
      },
    },
    {
      -logic_name => 'JobFactoryDefineMergedDataSet',
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::JobFactoryDefineMergedDataSet',
      -flow_into => {
	2 => 'DefineMergedDataSet'
      },
    },
    {
      -logic_name => 'JobFactoryPermissivePeakCalling',
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::JobFactoryPermissivePeakCalling',
      -flow_into => {
	100 => 'PermissiveSWEmbl'
      },
    },
    {
     -logic_name => 'MergeAlignments',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::MergeAlignments',
     -parameters => {
      	run_controls => 0,
     },
     -flow_into => {
# 	MAIN => 'RemoveDuplicateAlignments'
	MAIN => {
	  'RemoveDuplicateAlignments' => undef,
	  ':////accu?file_to_delete=[]' => { 
	    'file_to_delete' => '#bam_file_with_unmapped_reads_and_duplicates#'
	  }
	},
     },
     -rc_name => 'normal_monitored_2GB',
    },
    {
     -logic_name => 'RemoveDuplicateAlignments',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::RemoveDuplicateAlignments',
     -parameters => {
	run_controls => 0,
     },
     -flow_into => {
	  MAIN => [
	    'JobFactoryDefineMergedDataSet',
# 	    # Create bigwigs for replicates
# 	    WHEN('1' => { 'WriteBigWig' => INPUT_PLUS({ type => 'replicate' })})
	  ]
       },
#      -rc_name => '64GB_3cpu',
     -rc_name => 'normal_4GB_2cpu',
    },
    {
     -logic_name => 'MergeReplicateAlignments',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::MergeAlignments',
     -parameters => {
	run_controls => 0,
	permissive_peaks => $self->o('permissive_peaks')
      },
     -flow_into => {
#  	MAIN => 'RemoveDuplicateReplicateAlignments',
	MAIN => {
 	  'RemoveDuplicateReplicateAlignments' => undef,
	  ':////accu?file_to_delete=[]' => {
	    'file_to_delete' => '#bam_file_with_unmapped_reads_and_duplicates#'
	  }
	},
# 	# Create bigwigs for technical replicates
# 	MAIN => WHEN(
# 	  '1'  => { 'WriteBigWig' => INPUT_PLUS({ type => 'technical replicate' }) },
# 	),
     },
     -rc_name => 'normal_monitored_2GB',
    },
    {
     -logic_name => 'RemoveDuplicateReplicateAlignments',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::RemoveDuplicateAlignments',
     -parameters => {
	run_controls => 0,
     },
     -flow_into => {
	  MAIN => [
            'JobFactoryPermissivePeakCalling',
            # Create bigwigs for replicates
#             WHEN('1' => { 'WriteBigWig' => INPUT_PLUS({ type => 'replicate' })}),
	  ],
       },
#      -rc_name => '64GB_3cpu',
      -rc_name => 'normal_4GB_2cpu',
    },
    {
     -logic_name => 'WriteBigWig',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::RunWiggleTools',
     -rc_name    => 'normal_30GB_2cpu',
     -flow_into => {
	MEMLIMIT => 'WriteBigWig64GB',
      }
    },
    {
     -logic_name => 'WriteBigWig64GB',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::RunWiggleTools',
     -rc_name    => '64GB_3cpu',
    },
    {
      -logic_name    => 'PermissiveSWEmbl',
      -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunPeaks',
      -parameters => {
	peak_analysis => $self->o('permissive_peaks'),
      },
      -rc_name => 'normal_5GB_2cpu_monitored',
    },
    {
      -logic_name => 'CleanupFilesFromPermissiveSWEmblJobFan',
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::ErsaCleanup',
     -flow_into => {
	  MAIN => 'PreprocessIDR',
       },
    },
    {
     -logic_name => 'PreprocessIDR',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::PreprocessIDR',
     # They take an average of 7 seconds on the cttv dataset
     -batch_size => 30,
     -rc_name    => 'default',
     -parameters => {
	permissive_peaks => $self->o('permissive_peaks'),
      },
     },
   ];
}

1;

