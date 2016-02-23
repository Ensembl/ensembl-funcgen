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
	  '1' => [ 'TokenLimitedJobFactory' ],
	},
    },
    {   -logic_name => 'TokenLimitedJobFactory',
	-module     => 'Bio::EnsEMBL::Funcgen::Hive::TokenLimitedJobFactory',
	-meadow_type=> 'LOCAL',
	-flow_into => {
	  '2->A' => [ 'IdentifyAlignInputSubsets' ],
	  'A->1' => [ 'TokenLimitedJobFactory' ],
	},
    },
    {
      -logic_name => 'IdentifyAlignInputSubsets',
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::IdentifySetInputs',
      -meadow_type => 'LOCAL',#should always be uppercase
      -parameters => {
      set_type                     => 'InputSubset',
	feature_set_analysis_type    => 'peak',
	default_feature_set_analyses => $self->o('default_peak_analyses'),
	dataflow_param_names => ['no_idr'], 
      },
      -flow_into => {
	'2->A' => [ 'DefineResultSets'     ],
	'A->4' => [ 'CleanupCellLineFiles' ],
      },
      -analysis_capacity => 1, #For safety, and can only run 1 LOCAL?
      -rc_name => 'default',   #NA as LOCAL?
      # We don't care about these failing, as we expect them too
      -failed_job_tolerance => 100,
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
     #-module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
     #Inputs are either:
     #1 A set of controls and a series of InputSubset batches(replicates) which share the control
     #2 A series of unrelated InputSubset batches(replicates) which have no controls
    -parameters => {dataflow_param_names => ['no_idr']},#this is now batch flowed as we need it in DefineMergedOutputSet????
    -flow_into => {
      'A->2'  => ['PreprocessIDR'],#move to 2 'A->2' but make sure I flow this last for each control group
      '10'    => ['Preprocess_bwa_samse_control'],
      '11'    => ['Preprocess_bwa_samse_merged'],
      '12->A' => ['Preprocess_bwa_samse_replicate'],

      #Could add another alignment analysis here easily e.g.
      #'20'    => ['Preprocess_NEWALIGNER_control'],
      #'21'    => ['Preprocess_NEWALIGNER_merged'],
      #'22->A' => ['Preprocess_NEWALIGNER_replicate'],

      #Could support custom analyses like this
      #i.e. before we start the decimal blocks of
      #known resource optimised/branched analyses
      #Assuming we will never have more the 7 outflow analyses
      #allow_custom_branching would need to be an analysis level param
      #'3'    => ['Preprocess_custom_control'],
      #'4'    => ['Preprocess_custom_merged'],
      #'5->A' => ['Preprocess_custom_replicate'],

     },
     -analysis_capacity => 100,       #Will this just run 1 as LOCAL?
     -rc_name           => 'default', #NA as LOCAL?
    },

    ### PREPROCESS ANALYSES ###
    #Currently averaging ~22mins
    #But this depends entirely on the number and size of the reps in the control

    {
      -logic_name => 'Preprocess_bwa_samse_control',
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::PreprocessFastqs',
      -flow_into => {
	'2->A' => ['Run_bwa_samse_control_chunk'],
	'A->3' => ['MergeControlAlignments_and_QC'],
	},
      -batch_size => 1, #max parallelisation???
      -analysis_capacity => 200,
      -rc_name => '10gb_1cpu_staggered'
     },
     {
      -logic_name => 'Preprocess_bwa_samse_merged',
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::PreprocessFastqs',
      -parameters => {merge => 1},
      -flow_into =>
	{
	'2->A' => ['Run_bwa_samse_merged_chunk'],
	'A->3' => [ 'MergeAlignments_and_QC' ], #This is in Peaks.pm conf
	},
      -batch_size => 1, #max parallelisation???
      -analysis_capacity => 200,
      -rc_name => '10gb_1cpu_staggered'
     },
     {
      -logic_name => 'Preprocess_bwa_samse_replicate',
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::PreprocessFastqs',
      -flow_into => {
	'2->A' => ['Run_bwa_samse_replicate_chunk'],
	'A->3' => [ 'MergeReplicateAlignments_and_QC' ],
	},
      -batch_size => 1, #max parallelisation???
      -analysis_capacity => 200,
      -rc_name => '10gb_1cpu_staggered'
     },
    {
      -logic_name => 'Run_bwa_samse_control_chunk',
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::RunAligner',
      -batch_size => 20,
      #-analysis_capacity => 100,
      -rc_name => 'normal_10gb'
     },
    {
    -logic_name => 'Run_bwa_samse_merged_chunk',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::RunAligner',
     -batch_size => 20,
     #-analysis_capacity => 100,
     -rc_name => 'normal_10gb'
     },
    {
      -logic_name => 'Run_bwa_samse_replicate_chunk',
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::RunAligner',
      -batch_size => 20,
      #-analysis_capacity => 100,
      -rc_name => 'normal_10gb'
     },
    {
      -logic_name => 'MergeControlAlignments_and_QC',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::MergeQCAlignments',
     -parameters => {flow_mode => 'signal'},
     -flow_into => {
#        #2 is reserved for other Define DataSet flow (single vs multiple ResultSets).
#        # Although isn't really required as the branching is handled dynamically
#        #but let's keep it clean here for now.
#        'A->3'  => [ 'PreprocessIDR' ],
#        #alignment analyses encoded in blocks of 10
#        '10'    => ['Preprocess_bwa_samse_merged'],
#        '11->A' => ['Preprocess_bwa_samse_replicate'],
	  1 => 'JobFactorySignalProcessing',
       },
     -batch_size => 1, #max parallelisation
     -analysis_capacity => 200,
     -rc_name => '64GB_3cpu',
    },
    {
      -logic_name => 'JobFactorySignalProcessing',
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::JobFactorySignalProcessing',
      -flow_into => {
	'A->3'  => [ 'PreprocessIDR'                 ],
	'10'    => [ 'Preprocess_bwa_samse_merged'   ],
	'11->A' => [ 'Preprocess_bwa_samse_replicate'],
      },
      -meadow_type=> 'LOCAL',
    },
    {
      -logic_name => 'MergeAlignments_and_QC',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::MergeQCAlignments',
     -parameters => {flow_mode => 'merged'},
     -flow_into => { 2 => ['DefineMergedDataSet']},
     -batch_size => 1, #max parallelisation
     -analysis_capacity => 200,
     -rc_name => '64GB_3cpu',
    },
    {
     -logic_name => 'MergeReplicateAlignments_and_QC',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::MergeQCAlignments',
     -parameters => {
	flow_mode => 'replicate',
	permissive_peaks => $self->o('permissive_peaks')
      },
     -flow_into => {
	1 => 'JobFactoryPermissivePeakCalling'
     },
     -batch_size => 1, #max parallelisation
     -analysis_capacity => 200,
     -rc_name => '64GB_3cpu',
    },
    {
      -logic_name => 'JobFactoryPermissivePeakCalling',
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::JobFactoryPermissivePeakCalling',
      -flow_into => {
	'100' => 'run_SWEmbl_R0005_replicate'
      },
      -meadow_type=> 'LOCAL',
    },
    {
      -logic_name    => 'run_SWEmbl_R0005_replicate',  #SWEmbl permissive
      -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunPeaks',
      -parameters => {
	peak_analysis => $self->o('permissive_peaks'),
      },
      -analysis_capacity => 10,
      -rc_name => 'normal_5GB_2cpu_monitored',
    },
    {
     -logic_name => 'PreprocessIDR',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::PreprocessIDR',
     -batch_size => 30,
     -rc_name    => 'default',
     -parameters => {
	permissive_peaks => $self->o('permissive_peaks'),
      },
     },
   ];
}

1;

