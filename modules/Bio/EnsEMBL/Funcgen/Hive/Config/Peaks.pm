package Bio::EnsEMBL::Funcgen::Hive::Config::Peaks;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw );
use base ('Bio::EnsEMBL::Funcgen::Hive::Config::BaseDB');

sub pipeline_analyses {
  my $self = shift;

  return 
   [
    { #This is basically making sure the input file is sorted wrt genomic locations
	 -logic_name => 'PreprocessAlignments',
	 -module     => 'Bio::EnsEMBL::Funcgen::Hive::CollectionWriter',
	 -flow_into => 
      {
       #Not 3->A as we don't close these funnels
	   #These branch numbers must match the corresponding default config in DefineSets
       '3'   => [ 'run_SWEmbl_R015' ],      #Normal peaks
       '4'   => [ 'run_ccat_histone' ],     #CCAT
       '5'   => [ 'run_SWEmbl_R0025' ],     #Tight peaks 
       '6'   => [ 'run_SWEmbl_R0005_IDR' ], #Permissive peaks filtered by IDR max peaks value
       '7'   => [ ':////accu?file_to_delete_after_cell_line_has_been_processed=[file]' ], 
       '100' => [ 'run_peaks_custom' ],     #defined by feature_set_analysis
      },
	 -analysis_capacity => 100,
	 #Change this to hive_capacity as it may be competing with parallel peak jobs?
     -rc_name => 'default',
     #todo revise this to reserve tmp space as this is sorting the bed files
    },
    
 
  
  
    #we have an issue here as the fan changes from a slice fan, to an analysis fan
    #so branch 1 get's overloaded
    #do we need a factory here as we can't define fan on branch one and an accumulator
    #basically
    #should we just rebranch the collection confs to handle this 
    #or can we feed and wait for accumulators on the same branch?
    
    
      #Merge QC into this step? rather than flow to PeaksQC
      #Ww always want it run, but we don't really want to merge it into the PeakCaller
      #so merge it into RunPeaks?
      #We can do this for IDR as this runs on bed files i.e. before we do the write_output step
      #Until we do IDR (i.e. peak calling on separate replicates)
      #we will have to rely on peaks report
      #dowe do this for each set seprately or wait and report across all for comparison
  
      #what if only the PeaksQC failed, reload will enable us to rerun QC?
      #do we need added support for skip QC?
      #certainly need support for force_load if QC has failed
 
      #Do individual set Peak
 
    #Do not change these logic_names as they are used for dynamic dataflow/branching
    
 
	  {
	   -logic_name    => 'run_SWEmbl_R015',  #SWEmbl normal (histones/Dnase?)
	   -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunPeaks',
	   -analysis_capacity => 100,
	   -rc_name => 'normal_5GB_2cpu_monitored', # Better safe than sorry... size of datasets tends to increase...	   
	   
	  },

	  {
	   -logic_name    => 'run_ccat_histone', #CCAT
	   -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunPeaks', 
	   -parameters    =>
	    {process_file_types => ['significant_region'],
	     CCAT_parameters => {-chr_file => $self->o('data_root_dir').
	                           '/reference_files/CCAT/'.$self->o('species').'_'.$self->o('assembly').'.CCAT_chr_lengths.txt'}},
	   -analysis_capacity => 100,
	   -rc_name => 'normal_monitored_2GB', # CCAT does not need much?
	   #we were getting MEMLIMIT failures with default 200mb
	  
	  },
      {
	   -logic_name    => 'run_SWEmbl_R0025', #SWEmbl TF
	   -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunPeaks', 
	   -hive_capacity => 100,
	   -rc_name => 'normal_5GB_2cpu_monitored', # Better safe than sorry... size of datasets tends to increase...
	  # -flow_into => 
	  },
	  
	  
	  #Does this need to be shared with IDRPeaks
	  #actually making analysis names tied to branching
	  #means that we can't reuse the mechanism if we want to flow out differently
	  #Luckily this is not the case
	  
	  {
     -logic_name    => 'run_SWEmbl_R0005_IDR', #SWEmbl permissive IDR filtered
     -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunPeaks', 
     -parameters    => {'filter_max_peaks' => 1},
     -hive_capacity => 100,
     -rc_name => 'normal_5GB_2cpu_monitored', # Better safe than sorry... size of datasets tends to increase...
    # -flow_into => 
    },
    
    
    
	  #Add others in here e.g. macs etc
	  
      { 
       -logic_name    => 'run_peaks_custom', #Peaks caller not supported by config 
       -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunPeaks', # Peak caller defined by feature_set_analysis
       -hive_capacity => 100,
       -rc_name => 'normal_5GB_2cpu_monitored', # Better safe than sorry... size of datasets tends to increase...
       #-flow_into => 
      },
	 ];
}

1;
