
=head1 NAME

Bio::EnsEMBL::Funcgen::Hive::DefineInputSets

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::Hive::DefineResultSets;

use warnings;
use strict;
 
use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( scalars_to_objects );
use base ('Bio::EnsEMBL::Funcgen::Hive::BaseDB');

#Change this to DefineResultSets and drop InputSet!
#This will mean input_set.analysis moves to input_subset.
#input_set.replicate will move to result_set

#feature_set.input_set_id can be replaced with result_set_id
#this is only required for speed of experiment view and Source label zmenu rendering

#what are the implications for non seqencing input_sets?
#i.e. other flat file imports? e.g. segmentation
#This is non-critical and should really not be in the input_set/subset
#tables. Was just added there as this was a existing mechanism of import.

#should we rename input_subset as input_set

#what are the implications for tracking?

#what are the implications for exiting pipeline modules?
# IdentifySetInputs will have to change support for input_sets to result_sets

#Impact on existing import code. This is all based around input_set tracking.
#So this will need to be change to result_set. This needs to be overhauled for collections
#anyway, as this would currently store the bed file as an input_subset

#ResultSetAdaptor support

#How can we do this will minimal impact on web and release 74?


#Can't have no result_set mode now as this will break the link to the input_subset
#now input_set_id won't be in feature set.
#actually this is still fine for IDR sets, as these will never make it to dev.

sub fetch_input {   # fetch parameters...
  my $self = shift;
  $self->SUPER::fetch_input;
  $self->init_branching_by_analysis;
  
  my $isset_ids = $self->get_param_method('input_subset_ids',  'required');
  assert_ref($isset_ids, 'HASH', 'InputSubset dbIDs');
  
  $self->set_param_method('controls', delete $isset_ids->{controls});
  
  if(scalar(keys %$isset_ids) == 0){
    throw('No (signal) input_subset_ids defined. Must pass input_subset_ids hash of (optional) \'controls\''.
      ' and ResultSet name keys and input_subset_id arrayref values');
  }

  $self->get_param_method('alignment_analysis', 'required');
  #alignment_analysis is dataflowed explicitly from IdentifyAlignInputSubsets
  #to allow batch over-ride of the default/pipeline_wide value.
  #Could have put this in batch params, but it is only needed here
                    
  $self->get_param_method('merge_replicate_alignments', 'silent');

  #probably need to batch flow just -no_idr
  #as that will also be required in DefineMergedOutputSet? why?
  #Or can we just omit analysis to use default there?

  return;
}

sub _are_controls{
  my $ctrls        = shift; 
  my $all_controls = 1;
  
  foreach my $ctrl(@$ctrls){
    
    if(! $ctrl->is_control){
      $all_controls = 0;
      last;
    }
  }
  
  return $all_controls;
}

sub _are_signals{
  my $sigs     = shift; 
  my $all_sigs = 1;
  
  foreach my $sig(@$sigs){
    
    if($sig->is_control){
      $all_sigs = 0;
      last;
    }
  }
  
  return $all_sigs;
} 

sub run {   # Check parameters and do appropriate database/file operations... 
  my $self         = shift;
  my $helper       = $self->helper;
  my $rset_adaptor = $self->out_db->get_ResultSetAdaptor;
  my $issets       = $self->input_subset_ids;  
  my $ctrls        = [];  
  my $branch;
  
  #The branch will either be 2 for all, or may flip flop betwee 3 & for for idr and merged sets without controls
  
  my $align_lname = $self->alignment_analysis;
  my $align_anal  = &scalars_to_objects($self->out_db, 'Analysis',
                                                       'fetch_by_logic_name',
                                                       [$align_lname])->[0];
  
  #We don't use fetch_Set_input here as we are dealing with many InputSubsets
  
  if(@{$self->controls}){
    $ctrls = &scalars_to_objects($self->out_db, 'InputSubset',
                                                'fetch_by_dbID',
                                                $self->controls);
    if(! &_are_controls($ctrls)){
      throw("Found unexpected non-control InputSubsets specified as controls\n\t".
        join("\n\t", map($_->name, @$ctrls)));
    }
  }
                                             
  #Now we need to know whether to create a merged InputSet or a Replicate InputSet?
  #The set_names will already be set, but we actually need to set a rep value in the ResultSet
  #Just parse from the set_name, or set explicitly in the hash
  
  #validate all input_subsets before we start creating anything, such that rerunning the job will be clean
  
  #my %rsets = (2=> [], 3=> [], 4=>[]);#Withcontrols, Merged(no controls), IDR(no controls)
  #should this be done by init_branch_config?            
 
  my %rsets;
 
 
  #setname here might not actually be set name
  #$issets is actually:
  #{input_subset_ids => {$set_name => {$replicate => [dbID1]}}
  #                       controls => have been deleted above
  
  
  #why did we need these rep keys?
   
  foreach my $set_name(keys %$issets){
    my $parent_set_name = $set_name.'_'.$align_anal->logic_name;
    my $sigs = &scalars_to_objects($self->out_db, 'InputSubset',
                                                  'fetch_by_dbID',
                                                  $issets->{$set_name});
    if(! &_are_signals($sigs)){
      throw("Found unexpected controls in signal InputSubsets\n\t".
      join("\n\t", map($_->name, @$sigs)));
    }

    my $ftype        = $sigs->[0]->feature_type;
    my $is_idr_ftype = $self->is_idr_feature_type($ftype);
    my $merge_reps   = $self->merge_replicate_alignments;
    my $ctype        = $sigs->[0]->cell_type;
    my @rep_sets     = ($sigs);
    
    
    if($is_idr_ftype){
      
      if($merge_reps){
        #Merged control file should already be present
        #Merged signal file maybe present if GeneratePseudoReplicates has been run
        $branch = 'DefineMergedDataSet';
        
        
        #Merging of fastqs is normally done in PreprocessFastqs
        #But now we need to merge bams (or other)
        #Should we move this to PreprocessAlignments?
        #This normally expects the bams to be in merged if required.
        #It seem much more natural to do the merge here
        
        #Sanity check here the control file is available?
        #Using the ResultSet below?
        
        #Also do the merge handling below based on $merge_reps, and pulling back
        #the relevant rsets to get their alignment files.
        
        #Actually, probably best to do the merge in PostProcesIDR
        #as we have access to the rep ResultSets (and hence the relevant alignments files) there.
        
        
      }
      else{
        #RunIDR semaphore handled implicitly later 
        $branch = 'Preprocess_'.$align_lname.'_replicate';
        @rep_sets = map {[$_]} @$sigs;
      
        if(@$ctrls){
          $branch = 'Preprocess_'.$align_lname.'_control';
        }
      }
    }
    elsif(! @$ctrls){
      $branch = 'Preprocess_'.$align_lname.'_merged';
    }
    
    
    #So this is now creating the right sets, but how are we going to group the dbIDs/set_names
    
    
    #We need to add support in here for merging previously aligned replicate bams
    
    
    
    my $filter_format = $self->param_silent('bam_filtered') ? undef : 'bam'; 
    
    foreach my $rep_set(@rep_sets){ 
      my $rset_name = $parent_set_name.'_'.$align_anal->logic_name;
      my %rep_bams;
      
      if($is_idr_ftype){
             
        if($merge_reps){ #Iterate over all aligned sig ResultSets
          $rep_bams{rep_bams} = [];
          
          foreach my $rep(@$rep_set){
            #Pull back the rep Rset to validate and get the alignment
            #file for merging
            my $rep_rset_name = $rset_name.'_TR'.$rep->replicate;
            my $rset = $rset_adaptor->fetch_by_name($rep_rset_name);            
            #We could alternatively fetch them like this $rset_a->fetch_all_by_supporting_Sets.
            
            if(! defined $rset){
              $self->throw_no_retry("Could not find ResultSet for post-IDR merged ResultSet:\t".
                $rep_rset_name); 
            }
        
            #Validate controls are the same!
            #This should already have been done in PreprocessIDR
            #but probably a good idea to do here too
            #As we may get here by means other than PreprocessIDR?
            
            #Need to handle filtering here
            #at this point bam_filtered refers to the rep bams
            #not the output bam
            
            #Currently hardcoded for bam
            #although if we are merging can we/do we need to request sam?
            #what about clean up?
            push @{$rep_bams{rep_bams}}, 
              $self->get_alignment_file_by_ResultSet_formats($rset, ['bam'],
                                                             undef,  # control flag
                                                             undef,  # all_formats flag
                                                             $filter_format);
            
            
            
          }
        }
        else{ #only 1 in the $rep_set 
          $rset_name .= '_TR'.$rep_set->[0]->replicate;
        }
    
    
      #Run_BWA_and_QC_control will also need to know which way to branch
      #and the only way of doing that at present is implicit from the set name (TRN)
      #or of we do this above test again i.e. not_id grep broad_peak_feature_types
      #we could flow this explicitly? (although we need Run_BWA_and_QC_control to work independantly)
      #or just move that method to BaseSequenceAnalysis? (probably not worth it)
      #and logic/flow control will probably be different, dependant on context
      
      
      #Set the ref we want to push onto here?
      #Can't do this here as don't yet have anything to flow
      #but we do have the keys/grouping
      my $cache_ref;
  
  
      #change grouping here based on $run_reps and branch?
      #This will mean will mean we will have to break the dbIDs and set_names
      #structure, but this is fine  
  
  
      
      #if($branch =~ /_replicate$/o){ #Group IDR sets by merged/parent set
          
        
      if($run_reps){
        #branch can be control or replicate  
        
        #How can this be control? And isn't that merged?
        
        #This shoudl be data flown from IdentifyAlignInputSets
        $rsets{$branch}->{$parent_set_name} ||= [];  
        $cache_ref                            = $rsets{$branch}->{$parent_set_name};
      }
      else{
        #branch is only ever merged (no controls)
        $rsets{$branch}{merged} ||= [];
        $cache_ref                = $rsets{$branch}{merged};
      }
      

      #we don't create here as we want to validate the whole lot of input_subsets
      #before we create any
      #this will help with recovery

      push @$cache_ref,
       {-NAME          => $rset_name,
        #-FEATURE_CLASS        => 'result', # | dna_methylation',#This needs revising?
        #currently set dynamically in define_ResultSet
        -INPUT_SUBSETS => [@$rep_set, @$ctrls],
        -DBADAPTOR     => $self->out_db,
        -ANALYSIS      => $align_anal,
        -ROLLBACK      => $self->param('rollback'),#?
        -RECOVER       => $self->param('recover'), #?
        -CELL_TYPE     => $ctype,
        -FEATURE_TYPE  => $ftype,
        %rep_bams
       };
    }
  }
  
  #Now do the actual ResultSet generation and cache the output_ids on the correct branch
  my %batch_params = %{$self->batch_params};
  my %branch_sets;

  foreach my $branch(keys %rsets){
  
    foreach my $rset_group($rsets{$branch}){

      foreach my $rset(@$rset_group){
        
        my $rep_bams = delete($rset->{rep_bams});   
        $rset = $helper->define_ResultSet(%{$rset});
        
        if($rep_bams){
          #Check we don't already have the file from the Generate PseudoReps step
          #We will most likely already have filtered files by now, but we don't need to handle bam_filtered
          #as this will be done in PreprocessAlignments   
          my $merged_file = $self->get_alignment_file_prefix_by_ResultSet($rset).'.bam';
          
          if(! -f $merged_file || $self->param_silent('overwrite')){
            
            merge_bams($merged_file, $rep_bams, 
                       {remove_duplicates => 1,
                        sam_header        => $self->sam_header($rset->cell_type->gender)});
          }
        }
        
        if($branch =~ /_merged$/){# (no control) job will only ever have 1 rset    
          $self->branch_job_group($branch, [{%batch_params,
                                             dbID       => [$rset->dbID], 
                                             set_name   => [$rset->name],
                                             set_type    => 'ResultSet'}]);
        }
        else{
          #populate set_names/dbIDs groups for control and replicate branches
          $branch_sets{$branch}{$rset_group}{set_names} ||= [];
          $branch_sets{$branch}{$rset_group}{dbIDs}     ||= [];
          push @{$branch_sets{$branch}{$rset_group}{set_names}}, $rset->name;
          push @{$branch_sets{$branch}{$rset_group}{dbIDs}},     $rset->dbID;
  
          if($branch =~ /_replicate$/o){
            #Just create the individual fan output_ids
            $branch_sets{$branch}{$rset_group}{output_ids} ||= [];
            push @{$branch_sets{$branch}{$rset_group}{output_ids}}, {%batch_params,
                                                                     dbID     => $rset->dbID,
                                                                     set_name => $rset->name,
                                                                     set_type => 'ResultSet'};    
          }
        }
      }
    }
  }
  
  
  #Now flow job_groups of rsets to control job/replicate & IDR 
  foreach my $branch(keys %branch_sets){            
  
    if($branch =~ /_replicate$/){
      
      foreach my $rep_set(keys %{$branch_sets{$branch}}){ 
        #Add semaphore RunIDR job
        $self->branch_job_group($branch, $rep_set, 'RunIDR', 
                                [{%batch_params,
                                 dbIDs     => $branch_sets{$branch}{$rep_set}{dbIDs},
                                 set_names => $branch_sets{$branch}{$rep_set}{set_names},
                                 set_type  => 'ResultSet'}]);
      }           
    }       
    else{#($branch =~ /_control$/){    
      #Pick an arbitrary set for access to the controls 
      my ($any_group) = keys(%{$branch_sets{$branch}});
   
      $self->branch_job_group($branch, 
                              [{%batch_params,
                               result_set_groups => $branch_sets{$branch},
                               set_type  => 'ResultSet',
                               dbID      => $branch_sets{$branch}{$any_group}{dbIDs}->[0],
                               set_name  => $branch_sets{$branch}{$any_group}{set_names}->[0]}]);   
    }
  }
  
  return;
}



sub write_output {  # Create the relevant jobs
  shift->dataflow_job_groups;
  return;
}

1;