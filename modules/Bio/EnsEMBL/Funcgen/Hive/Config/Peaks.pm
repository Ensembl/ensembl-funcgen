package Bio::EnsEMBL::Funcgen::Hive::Config::Peaks;

use strict;
use warnings;
use base ('Bio::EnsEMBL::Funcgen::Hive::Config::Base');

sub pipeline_analyses {
  my $self = shift;

  return 
   [
      { #This is basically making sure the input file is sorted wrt genomic locations
	    -logic_name => 'PreprocessAlignments',
	    -module     => 'Bio::EnsEMBL::Funcgen::Hive::CollectionWriter',
	    -flow_into => {
	  #Not 3->A as we don't close these funnels
	      #These branch numbers must match the corresponding default config in DefineSets
	  '3'   => [  'run_SWEmbl_R015' ],      #Normal peaks
	  '4'   => [  'run_ccat_histone' ],     #CCAT
	  '5'   => [  'run_SWEmbl_R0025' ],     #Tight peaks 
	  '6'   => [  'run_SWEmbl_R0005_IDR' ], #Permissive peaks filtered by IDR max peaks value
	  '7'   => [ ':////accu?file_to_delete_after_cell_line_has_been_processed=[file]' ], 
	  '100' => [  'run_peaks_custom' ],     #defined by feature_set_analysis
	},
      },
      #Do not change these logic_names as they are used for dynamic dataflow/branching
      {
	-logic_name    => 'run_SWEmbl_R015',  #SWEmbl normal (histones/Dnase?)
	-module        => 'Bio::EnsEMBL::Funcgen::Hive::RunPeaks',
	-rc_name => 'normal_5GB_2cpu_monitored',
      },
      {
	-logic_name    => 'run_ccat_histone', #CCAT
	-module        => 'Bio::EnsEMBL::Funcgen::Hive::RunPeaks', 
	-parameters    => {
	  process_file_types => ['significant_region'],
	  CCAT_parameters => {
	    -chr_file => $self->o('chromosome_file')
	  }
	},
	-rc_name => 'normal_monitored_4GB',
      },
      {
	-logic_name    => 'run_SWEmbl_R0025', #SWEmbl TF
	-module        => 'Bio::EnsEMBL::Funcgen::Hive::RunPeaks', 
	-rc_name => 'normal_5GB_2cpu_monitored',
      },
      {
	-logic_name    => 'run_SWEmbl_R0005_IDR', #SWEmbl permissive IDR filtered
	-module        => 'Bio::EnsEMBL::Funcgen::Hive::RunPeaks', 
	-parameters    => {'filter_max_peaks' => 1},
	-rc_name => 'normal_5GB_2cpu_monitored',
      },
      # Add others in here e.g. macs etc
      { 
       -logic_name    => 'run_peaks_custom', #Peaks caller not supported by config 
       -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunPeaks', # Peak caller defined by feature_set_analysis
       -rc_name => 'normal_5GB_2cpu_monitored',
      },
   ];
}

1;
