
=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::Hive::Config::DefineMergedDataSets;

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 CONTACT

    Please contact ensembl-dev@ebi.ac.uk mailing list with questions/suggestions.

=cut


package Bio::EnsEMBL::Funcgen::Hive::Config::DefineMergedDataSets;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Funcgen::Hive::Config::BaseSequenceAnalysis');
# All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly

#todo this change to Signal/Wiggle.pm config 
#with DefineSetInputs and DefineSets as shared analyses

#We really want to be able to genericise this such that we can swap out
#any given Runnable/DB. Move all specific options to parameters sections
#specific for that analysis


=head2 default_options

    Description : Implements default_options() interface method of
    Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that is used to initialize default options.

=cut


sub default_options {
  my ($self) = @_;
 
  #Define any optional params here as undef
  #If they are mandatory for one analysis, but optional for another
  #will have to catch that in the runnable
  
 
  return 
    {
     %{$self->SUPER::default_options},        
     run_DefineMergedDataSet => 1, 
  };
}


=head2 pipeline_wide_parameters

    Description : Interface method that should return a hash of pipeline_wide_parameter_name->pipeline_wide_parameter_value pairs.
                  The value doesn't have to be a scalar, can be any Perl structure now (will be stringified and de-stringified automagically).
                  Please see existing PipeConfig modules for examples.

=cut


#Can we move some of these to Base.pm config?

#sub pipeline_wide_parameters {
#  my $self = shift;
  
   #Deal with batch_params first as this may not have been subtituted yet
   #and will evaluate to a string:
   #Can't use string ("#:subst batch_params:#") as an ARRAY ref while "strict refs"
     
  # my $batch_params = 
  #    [ ( ref($self->o('batch_params')) ? 
  #        @{$self->o('batch_params')} : () ),
  #      #These allow override of defaults  
  #      'result_set_analysis',
  #      'feature_set_analysis',
  #    ];
             
#  return 
#    {
#     %{$self->SUPER::pipeline_wide_parameters}, 
#    };
#}


=head2 pipeline_create_commands

    Description : Implements pipeline_create_commands() interface method of
      Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that lists the commands
      that will create and set up the Hive database.

=cut

## WARNING!!
## Currently init_pipeline.pl doesn't run this method when a pipeline is created with the -analysis_topup option



=head2 pipeline_analyses

    Description : Implements pipeline_analyses() interface method of
      Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that defines the structure of the pipeline: analyses, jobs, rules, etc.


=cut



#(2) The hive_capacity worker limiting mechanism that was in place for years is going to change slightly.
#The old meanings of the values were:
#	negative value  : checking is switched off in this particular analysis
#	zero		: no workers will be allowed to take this analysis
#	missing value   : sets the default, which is 1
#	positive value	: (including the 1 by default) this analysis is limited by the given value
#This was counter-intuitive, because by default there was already a limit, which people had to raise or switch off if they needed.
#Now, since we also have an alternative mechanism (analysis_capacity), I'd like to make both mechanisms "off" by default
#(missing value will mean "checking is switched off").
#
#So please, please, please - check whether your pipeline *RELIES* on the current default behaviour
#(which is "if you don't specify -hive_capacity, the system will set it to 1
#and will not allow any workers of other analyses to take any jobs while this analysis is running").
#
#If you know/find out that your pipeline *RELIES* on this functionality, please explicitly set -hive_capacity => 1 in the corresponding analyses.
#This will make your pipeline compatible with future releases of the Hive system.



sub pipeline_analyses {
  my $self = shift;

  return [
    {
     -logic_name => 'IdentifyMergedResultSets',
	 -module     => 'Bio::EnsEMBL::Funcgen::Hive::IdentifySetInputs',	  
	 -meadow_type => 'LOCAL',#should always be uppercase
	  
	 #general parameters to pass to all jobs, use_tracking_db?
	 -parameters => {set_type        => 'result_set',
	                 #These are similar to batch_params but only flow
	                 #to the next analysis
	                 },
	                 #dataflow_param_names => ['recover']}, 
	                 #this is redundant as it is a batch_param_name!

	  	  
	  	  
  	  #Do we need to change the param names here?
  	  #When we merge confs, this will mean that
  	  #subsequent IdentifySetInputs job will get all the params
  	  #from the original init?
  	  #are the input_ids ignored when we merge? omitted
	  	  
  	 #These are always created, even if we have used analysis_top_up
  	 #hence at the start of the next conf, we will get input_ids which are data flown
  	 #Do we even want to redo this?
  	 # as this forces use to funnel
  	 #when what we really want is to data flow directly from batch job to batch job
  	 #Do we need a sparse shared batch job analysis?
  	 
  	 #regardless, the way to fix this is probably to remove these from here
  	 #and use seed_pipeline after init_pipeline
  	 #seed pipeline func could automatically print out the jobs created? 
	 #can seed pipeline use the existing meta table params?
	 #So init will do the validation and seed pipeline will simply
	 
	 #This is actually fine here, as this is the very first 
	 #shared IdentifySetInputs analysis
	 #Although we should probably make this consistent across the pipelines
	 #How do we catch mistyped options?
	 #Need to put these in pipeline_wide_parameters otherwise they won't be caught?
	 #This won't give us the ability to translate between input_sets and set_names
	 #This won't matter as the specificity will be implied by the logic_name
	 #although not in the meta table
	 #Not too happy about separating the input_id definition from this config
	 
	 #Putting them in pipelinewide params offers not advantage
	 #and may prevent using seed_pipeline again, as it would require 
	 #init_pipeline to be run again, probably unsafe for various reasons
	 #including over-writing of meta params
	 #let's write our own seed_pipeline script
	 #which calls hive seed_pipeline, but we catch opts better
	 #and optionally print input_id as output(all or newly seeded)
	 #this would also handle array args!
	 
	 
	  	  
	 #-input_ids => #Single input_id containing params to identify InputSets
	 #  [{
	 #    'cell_types'    => $self->o('cell_types'), 
     #    'feature_types' => $self->o('feature_types'), 
     #    'set_names'     => $self->o('input_sets'),#This is mutally exclusive to others
     #    'set_ids'       => $self->o('input_set_ids'),
     #    'experimental_groups' => $self->o('experimental_groups'),
     #    'states'              => $self->o('states'),
	 #	   #, 'result_file' => $self->o('result_file')  
	 #   }],
			 	 
	 -flow_into => {		 
      '2' => [ 'DefineMergedDataSet' ],
      #Not 2->A as we don't close the funnel for this pipeline
      #'A->1' => [ ],#Nothing else waiting for this job
      #Use branch 2 here so we can flow input_ids
      #down branch 1 back bone later to a merge step
      #although we will need to explicitly data flow
      #Is there a post MergeCollections step?
      
      
      
      #IDR/Replicate Implementation
      #We will need additional semaphored data_flow to run_peak_replicates here
      #to support replicate calling and IDR QC
      #This needs to wait for all replicate sets to be peak called
      #but ignore the rest(we don't want to wait for merged replicate peaks)
      #there is no way of doing this without splitting the whole analysis tree from
      #here i.e. DefineOutputSet?!!!
      #it's clear that implementing IDR optionally along side merged peak calling 
      #in the same pipeline is going to be tricky
      #This would need to flow to another link analysis(if >1 is supported)
      #Something like Run_IDR_QC
      #That link analysis would pick up the groups of input_set_ids
      #do the IDR across the linked FeatureSets and the submit another merged 
      #peak job.
      #This would mean duplicating the peak analyses as we can't loop back due to semaphore
      #IdentifyFeatureSets(or more likely PeaksReport) would also have to be semaphored from here
      #and receive the output of all peak jobs from both originally merged and IDR jobs?
      #No, we can handle this with status entries
      #Further more PreprocessALignment should not data flow to WriteCollections if it is 
      #a replicate
      
      #How can we semaphore an analysis(PeaksReport) which doesn't exist in this config?
      #would have to do this through the link analysis
      #This would wait for everything hanging of DefineOutputSet before commencing 
      #the link analysis, which would then run the IDR and merged peaks, before doing the final 
      #PeaksReport
      #This would also wait for Collections! So we would have to detach that analysis somehow
      #or split DefineOutputSet into DefineReplicateOutputSet and DefineOutputSet
      #VERY COMPLICATED!!!
      
      },
	
	   
	
			 
	  -analysis_capacity => 10,
      -rc_name => 'default',
    },
	  
    {
      #This will define a ResultSet usign the InputSet from above
      #or if used with -analysis_topup and the Peaks.pm config
      #will define a FeatureSet and DataSet with a ResultSet as support
      #along with the InputSet from above 
      
      #How will this know to make just a ResultSet or the full DataSet?
      
      
     -logic_name => 'DefineMergedDataSet', 
	   -module     => 'Bio::EnsEMBL::Funcgen::Hive::DefineDataSet',
	   -parameters => 
	    {#Need feature here, as we DefineDataSetmay want some for result_set too?
	     #although this is now currently impossible since we are removing InputSet
	     #but let's keep this generic for now
	     default_feature_set_analyses => $self->o('default_peak_analyses'),
	     feature_set_analysis_type    => 'peak',
	     #not $self->o('peak_analysis') as we want to batch flow this
	     #set_type=result_set could go here, but this is entirely dependant on what
	     #flows into DefineMergedDataSet, so leave that to IdentifyResultSets
	     #although this will likely ahve already been created
	     #rollback                     => $self->o('rollback'), 
	    },
	  #Add in default params from above!
	  
	  
	  
	  #-input_ids => [],#These will be explicitly flowed from DefineSetInputs
	 #Would be nice to use names here instead of dbIDs, but we can't guarantee
	 #this until the result set nameing convention is changed, and we make name a unique field
			 	 
	   -flow_into => 
	     {
		  '1' => [ 'PreprocessAlignments' ],
		   #Not 1->A as we don't close the funnel for this part of the pipeline
		 },
		 
	   -analysis_capacity => 10,
       -rc_name => 'default',
	   #this really need revising as this is sorting the bed files
	   #Need to change resource to reserve tmp space
	   
	   #Not having a -failed_job_tolerance here is causing the beekeeper to 
	   #exit, especially as there is no -max_retry_count set either
	   
    },
    
    
    #This needs changing to take a result_set_id and not create the sets
    #shall we merge this into the next step
    #No we still need this as this batches the jobs further
    #This will simply
    
     {
	   -logic_name => 'PreprocessAlignments',#was 'SetUp',
	   #This is basically making sure the input file is sorted wrt genomic locations
	   -module     => 'Bio::EnsEMBL::Funcgen::Hive::CollectionWriter',
	   
	   #Need to maintain this here as will not be updated by -analysis_topup
	   -parameters => 
	    {
	     feature_formats         => ['bam', 'bed'],
	     
	     #Would be nice to do something like this in each of the topup configs
	     #However, topped up analyses don't have their parameter updated!
         #-parameters => { feature_formats => '#expr([ @{#feature_formats#}, "bam"])expr#'},
	     
	     
	     #Would really like to build this array via the topped up confs
	     #i.e. each conf will add it's own formats
	     #$self->o doesn't allow inheritance style building of arrays
	     #let's try it with # variable # notation?     
	     peak_branches           => $self->o('peak_branches'),
	    },
		
	  	
		
	   #No -flow_into defined as this will be handled by -analysis_topup
	   # See Peaks and Collections conf.			 	 
			 
	   -analysis_capacity => 10,
       -rc_name => 'normal_2GB',
	   #this really need revising as this is sorting the bed files
	   #Need to change resource to reserve tmp space
	  },
    
	 
  ];
}




1;
