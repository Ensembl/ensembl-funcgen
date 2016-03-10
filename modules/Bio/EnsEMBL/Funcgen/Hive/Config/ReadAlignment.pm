package Bio::EnsEMBL::Funcgen::Hive::Config::ReadAlignment;

use strict;
use warnings;
use Data::Dumper;
use base ('Bio::EnsEMBL::Funcgen::Hive::Config::BaseSequenceAnalysis');

sub default_options {
  my $self = shift;
  return {
    %{$self->SUPER::default_options},

      #Size of each sequence chunk to be aligned (nbr of reads * 4)
      #fastq_chunk_size      => 16000000, #This should run in 30min-1h
      fastq_chunk_size      =>   1000000,
      alignment_analysis    => 'bwa_samse',
      bwa_samse_param_methods     => ['sam_ref_fai'],
      fastq_root_dir => $self->o('data_root_dir').'/fastq',
   };
}

sub pipeline_wide_parameters {
    my $self = shift;
    return {
      %{$self->SUPER::pipeline_wide_parameters},

      #Size of each sequence chunk to be aligned (nbr of reads * 4)
      fastq_chunk_size      => $self->o('fastq_chunk_size'),   #Change to batch specific
      alignment_analysis    => $self->o('alignment_analysis'), #Nope we may want this to be batch specific!
      aligner_param_methods => $self->o('bwa_samse_param_methods'),
      #This is stricly not required anymore as we use the local_url from the tracking tables
      fastq_root_dir      => $self->o('fastq_root_dir'),
      #This will should be set to one in downstream config
      can_PreprocessIDR              => 0,
      can_run_SWEmbl_R0005_replicate => 0,
      can_DefineMergedDataSet       => 0,
    };
}

sub pipeline_analyses {
  my $self = shift;
  return
   [
    @{$self->SUPER::pipeline_analyses}, #To pick up BaseSequenceAnalysis-DefineMergedOutputSet

    {   -logic_name => 'JobPool',
	-module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
	-wait_for    => 'JobPool',
	-flow_into => {
	  MAIN => 'TokenLimitedJobFactory',
	},
    },
    {   -logic_name => 'TokenLimitedJobFactory',
	-module     => 'Bio::EnsEMBL::Funcgen::Hive::TokenLimitedJobFactory',
	-meadow_type=> 'LOCAL',
	-flow_into => {
	  '2->A' => 'IdentifyAlignInputSubsets',
	  'A->1' => 'TokenLimitedJobFactory',
	},
    },
    {
      -logic_name => 'IdentifyAlignInputSubsets',
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::IdentifySetInputs',
      -meadow_type => 'LOCAL',
      -parameters => {
      set_type                     => 'InputSubset',
	feature_set_analysis_type    => 'peak',
	default_feature_set_analyses => $self->o('default_peak_analyses'),
	dataflow_param_names => ['no_idr'], 
      },
      -flow_into => {
	'2->A' => 'DefineResultSets',
	'A->4' => 'CleanupCellLineFiles',
      },
    },
    {
      -logic_name => 'CleanupCellLineFiles',
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::ErsaCleanup',
      -meadow_type=> 'LOCAL',
    },
    {
     -logic_name => 'DefineResultSets',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::DefineResultSets',
     -meadow     => 'LOCAL',
    -flow_into => {
      '2->A' => 'FixResultSetsExperimentIds',
      'A->3' => 'Preprocess_bwa_samse_control',
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
	sql => qq(
	  update result_set, cell_type, experiment 
	  set result_set.experiment_id = experiment.experiment_id, result_set.replicate = #replicate#
	  where cell_type.cell_type_id=experiment.cell_type_id and cell_type.cell_type_id=result_set.cell_type_id and experiment.feature_type_id=result_set.feature_type_id
	  and result_set_id = #dbID#
	  ),
	  db_conn => '#out_db_url#'
	},
      -meadow     => 'LOCAL',
      -analysis_capacity => 1,
      -batch_size => 1000,
    },
    {
      -logic_name => 'Preprocess_bwa_samse_control',
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::PreprocessFastqs',
      -flow_into => {
	'2->A' => 'Run_bwa_samse_control_chunk',
	'A->3' => 'MergeControlAlignments',
	},
      -rc_name => '10gb_1cpu'
     },
     {
      -logic_name => 'Preprocess_bwa_samse_merged',
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::PreprocessFastqs',
      -parameters => {merge => 1},
      -flow_into =>
	{
	'2->A' => 'Run_bwa_samse_merged_chunk',
	'A->3' =>  'MergeAlignments',
	},
      -rc_name => '10gb_1cpu'
     },
     {
      -logic_name => 'Preprocess_bwa_samse_replicate',
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::PreprocessFastqs',
      -flow_into => {
	'2->A' => 'Run_bwa_samse_replicate_chunk',
	'A->3' => 'MergeReplicateAlignments' ,
	},
      -rc_name => '10gb_1cpu'
     },
    {
      -logic_name => 'Run_bwa_samse_control_chunk',
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::RunAligner',
      -rc_name => 'normal_10gb'
     },
    {
    -logic_name => 'Run_bwa_samse_merged_chunk',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::RunAligner',
     -rc_name => 'normal_10gb'
     },
    {
      -logic_name => 'Run_bwa_samse_replicate_chunk',
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
	  'MAIN' => 'JobFactorySignalProcessing',
       },
     -rc_name => '64GB_3cpu',
    },
    {
      -logic_name => 'JobFactorySignalProcessing',
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::JobFactorySignalProcessing',
      -flow_into => {
	'A->3'  => 'PreprocessIDR',
	'10'    => 'Preprocess_bwa_samse_merged' ,
	'11->A' => 'Preprocess_bwa_samse_replicate',
      },
      -meadow_type=> 'LOCAL',
    },
    {
      -logic_name => 'JobFactoryDefineMergedDataSet',
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::JobFactoryDefineMergedDataSet',
      -flow_into => {
	2 => [ 'DefineMergedDataSet' ]
      },
      -meadow_type=> 'LOCAL',
    },
    {
      -logic_name => 'JobFactoryPermissivePeakCalling',
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::JobFactoryPermissivePeakCalling',
      -flow_into => {
	'100' => [ 'run_SWEmbl_R0005_replicate' ]
      },
      -meadow_type=> 'LOCAL',
    },
    {
     -logic_name => 'MergeAlignments',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::MergeAlignments',
     -parameters => {
#       flow_mode => 'merged',
      	run_controls => 0,
     },
     -flow_into => {
	MAIN => 'JobFactoryDefineMergedDataSet'
     },
     -rc_name => '64GB_3cpu',
    },
    {
     -logic_name => 'MergeReplicateAlignments',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::MergeAlignments',
     -parameters => {
	run_controls => 0,
	permissive_peaks => $self->o('permissive_peaks')
      },
     -flow_into => {
	MAIN => 'JobFactoryPermissivePeakCalling'
     },
     -rc_name => '64GB_3cpu',
    },
    {
      -logic_name    => 'run_SWEmbl_R0005_replicate',  #SWEmbl permissive
      -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunPeaks',
      -parameters => {
	peak_analysis => $self->o('permissive_peaks'),
      },
      -rc_name => 'normal_5GB_2cpu_monitored',
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

