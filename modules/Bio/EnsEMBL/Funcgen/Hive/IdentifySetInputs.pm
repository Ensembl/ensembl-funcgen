
=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.


=head1 NAME

Bio::EnsEMBL::Funcgen::Hive::IdenetifySetInputs;

=head1 DESCRIPTION

This module simply takes a list of paramters and a 'set_type' to identify 
Sets used as inputs for various parts of the analysis pipeline.

=cut

package Bio::EnsEMBL::Funcgen::Hive::IdentifySetInputs;

use warnings;
use strict;

use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( scalars_to_objects );
use base ('Bio::EnsEMBL::Funcgen::Hive::BaseDB');



#todo -slice_import_status?


my %set_adaptor_methods = 
(
  InputSubset => 'get_InputSubsetAdaptor',
  InputSet    => 'get_InputSetAdaptor',
  ResultSet   => 'get_ResultSetAdaptor',
  FeatureSet  => 'get_FeatureSetAdaptor',
  #data_set?
);


#Accessors to catch compile time errors and avoid uncaught typos

#These are now all injected by set_param_method below
#sub sets {                return $_[0]->param('sets');                }
#sub set_type {            return $_[0]->param_required('set_type');   }
#sub set_adaptor {         return $_[0]->param('set_adaptor');         }
#sub constraints_hash {    return $_[0]->param('constraints_hash');    }

#These are now injected by get_param_method via process_params below
#sub feature_types {       return $_[0]->param('feature_types');       }
#sub cell_types {          return $_[0]->param('cell_types');          }
#sub experimental_groups { return $_[0]->param('experimental_groups'); }
#sub analyses {            return $_[0]->param('analyses');            }
#sub set_ids {             return $_[0]->param('set_ids');             }
#sub set_names {           return $_[0]->param('set_names');           }
#sub states {              return $_[0]->param('states');              } 

#TODO 
#1 Check we have experiment/experimetal_group support for the composable queries
#  validate set type if not all SetAdaptors support this.
#2 Could do with a way of listing embargoed, when no_write is defined
#3 Implement is_InputSubset_downloaded in run
#4 Genericse no_idr into something like merge_reps_only? #This is aimed at restricting to non-idr sets
# and is not yet implemented. We really need a difference flag here. As no_idr is for turning off idr
# not set selection and maintaining idr. Change this to run_idr as a boolean? This is still a bit ambiguous
#5 review force/ignore_embargoed
#6 Review throw string building once we support download


sub fetch_input {   # fetch parameters...
  my $self = shift;
  $self->SUPER::fetch_input; 
  #$self->helper->debug(3, "IdentifySetInput params:\t", $self->input_job->{_param_hash});
    
  my $set_type = $self->get_param_method('set_type', 'silent'); 
  
  #This auto detection prevents having to do this in the SeedPipeline environment func
  if(! defined $set_type){
    my $lname = $self->analysis->logic_name;  
    
    foreach my $stype( qw( ResultSet InputSubset ) ){
      #no input_set support here    
      #(my $lname_type = $stype) =~ s/_//o;
    
      if($lname =~ /$stype/){
        
        if($set_type){ #Throw if we match is not unique
          throw("Cannot auto detect unique set_type from analysis logic name:\t$lname\n",  
                'Please defined the set_type analysis_parameter or rename analysis.');
        }
        
        $set_type = $stype;
      }
    }
    
    if(! defined $set_type){
      throw("Cannot auto detect unique set_type from analysis logic name:\t$lname\n",  
            'Please defined the set_type analysis_parameter or rename analysis.'); 
    } 
  }
    
  #Need to define the methods here, as we call them for all set types in run
  my $force_embargoed  = $self->get_param_method('force_embargoed',  'silent');
  my $ignore_embargoed = $self->get_param_method('ignore_embargoed', 'silent'); 
  #This currently still supports input_sets
 
 
  #set types in config are table case, where as the are some class case references in here
  #didn't we decide on class case? Yes!
 
  if($set_type =~ /^Input(Sub)?[sS]et$/) {
    $self->set_param_method('validate_InputSubset_tracking', 1); #for convenience/speed in loop below
        
    if($force_embargoed && $ignore_embargoed){
      throw('force_embargoed and ignore_embargoed are mutually exclusive parameters. '.
        'Please omit one or both');     
    }
    
    #Get feature_set analysis for control validation in run
    
  
    #but we should also be mindful of a peak_analysis over-ride?
    #but need to use generic params here?

    #we could replace this with control_mandatory_ftypes
    
    #Genericise feature_set_analysis_type to feature_set_analysis
    my $anal_type = $self->param_required('feature_set_analysis_type');     
    $self->set_param_method('feature_set_analysis', 
                            $self->param_silent($anal_type.'_analysis'));
    
    if(! defined $self->feature_set_analysis){
      $self->get_param_method('default_feature_set_analyses', 'silent');     
 
      if(! defined $self->default_feature_set_analyses){
        throw("Please define -${anal_type}_analysis or add to default_feature_set_analyses".
          " in the default_options config");
      }    
    }
  }
  
  if(! exists $set_adaptor_methods{$set_type}){
    throw("The -set_type $set_type is not supported by IdentifySetInputs.\n".
          "Valid options are ".join(' ', keys(%set_adaptor_methods)) );      
  }
  else{
    #Set the method name value to reset the value as the actual adaptor
    my $method = $set_adaptor_methods{$set_type};
    $self->set_param_method('set_adaptor', $self->out_db->$method, 'required');
  }
  
  
  $self->process_params([qw(set_names set_ids only_replicates allow_no_controls 
                            identify_controls control_experiments)], 1);#optional
              
  if($set_type eq 'InputSubset'){                      
    if($self->identify_controls && $self->allow_no_controls){
      throw("Mutually exclusive parameters have been specified:\t".
        "identify_controls & allow_no_controls\nPlease omit one of these.");  
    }
    
    if($self->control_experiments  && $self->allow_no_controls){
      throw("Mutually exclusive parameters have been specified:\t".
        "control_experiments & allow_no_controls\nPlease omit one of these.");  
    }
    
                    
    #if($self->control_experiments  && $self->identify_controls){
    #  throw("Mutually exclusive parameters have been specified:\t".
    #    "control_experiments & identify_controls\nPlease omit one of these.");  
    #}
    
    #if(! ($self->control_experiments || 
    #      $self->identify_controls ||
    #      $self->allow_no_controls) ){
    #  throw("For set_type input_subset, one of the following must be defined:\t".
#        'control_experiments identify_controls or allow_no_controls');
#    }
#Not strictly true we may only have subsets which are associated with experiments which have controls
  }
  elsif($self->control_experiments || 
        $self->identify_controls ||
        $self->allow_no_controls){
    throw('You have specified a parameter which is not appropriate for '.$set_type.
      "s:\n\tcontrol_experiment or identify_controls or allow_no_controls");
  }  
  
  
  #Set an internal handle_replicates param
  #so we don't try and do this for non input_sub/sets.
  my $only_reps = 0;
  
  if($self->only_replicates){
    #Filter replicates from generic fetch_all query
    #and validate specific id/name queries return only reps
  
    #No this is wrong, this should only be for ResultSets
    #Or can we support both for now
    #InputSets will go away
  
  
    if($set_type eq 'InputSet'){
      $only_reps = 1;    
    }
    else{
      throw("only_replicates param is only valid for InputSets, not ${set_type}s.\n".
        'Please omit only_replicates config or set the set_type to InputSet');
    }
  }
  
  #Revise only_replicates/handle_replicates
  #Currently this is ambiguous
  #why did this get renamed to handle_replicates for the run step?
  #currently only used in IdentifyReplicateResultSets
  #but was originally aimed at handling InputSets too, although these will now be dropped
  #Current usage is also a proxy for flowing permissive_peaks to 
  #run_SWEmbl_R0005_replicate
  #Although why do we even need this, can't RunPeaks pick up run_idr and is_idr_feature_type?
  #RunPeaks, should have to know about IDR!
  
  $self->set_param_method('handle_replicates', $only_reps);
  
  #Parse comma separated lists of filters into arrayrefs of string or objects
  $self->set_param_method('constraints_hash',
                          $self->process_params([qw(feature_types cell_types states
                                                    analyses experiments experimental_groups)],
                                                1, 1)  #optional/as array flags
                          );      

  #Catch mutally exclusive filter params  
  
  if($self->feature_types ||
     $self->cell_types ||
     $self->experiments ||
     $self->experimental_groups ||
     $self->analyses ||
     $self->states){ 
     #all these are OR filters except states which is an AND filter
      
    if($self->set_names || $self->set_ids){
      throw('You have specified mutually exclusive filter params for the '.
            "IdentifySetInputs analysis\nPlease specify restrict to ".
            '-set_name or -set_ids or a combination other filters '.
            '(e.g. -experimental_groups -experiments -feature_types -cell_types -analyses -states');
    }
  }
  elsif(! ($self->set_names || $self->set_ids)){
 
    throw('You must specifiy some IdentifySetInputs filter params either '.
          '-input_sets or -set_ids or a combination of '.
          '-feature_types -cell_types -experiments -experimental_groups -states -analyses');
  }
  elsif($self->set_names && $self->set_ids){
    throw('You have specified mutually exclusive filter params for the '.
            "IdentifySetInputs analysis\nPlease specify restrict to ".
            'set_names or -set_ids or a combination other filters '.
            '(e.g. -experimental_groups -experiments -feature_types -cell_types -analyses -states');  
  }


  #$self->init_branch_config(1, 1);Ê#No longer needed
  #Flags are:
  #Optional branch config, as we only need this for IdentifyAlignInputsets 
  #Validate branch_key_method if we do have config

  return;
}



sub run {   # Check parameters and do appropriate database/file operations... 
  my $self = shift;
  #Grab dataflow_params here as these params will not 
  #change as we iterate through sets below
  #These will only be flowed to the next job
  #where as batch_params flow across all jobs in this seeded batch
  #(for those that support/require batch params)
  my $dataflow_params = $self->dataflow_params(1);#optional flag
  my $batch_params    = $self->batch_params; 
  my $handle_reps     = $self->handle_replicates;
  my $set_adaptor     = $self->set_adaptor;
  my $set_type        = $self->set_type;
  my $no_write        = $self->param_silent('no_write'); 
  my $sets  = [];
  my $throw = '';
  my (@failed);
  
  ### FETCH SETS ###
 
 
  #Todo
  #1 Changed @failed to build $throw directly with more contextual info
  
  #For set_ids and set_names, catch undef return types
  if($self->set_ids || $self->set_names){
    my ($ids_or_names, $method);
    
    if($self->set_ids){
      $method       = 'fetch_by_dbID'; 
      $ids_or_names = $self->set_ids; 
    }
    else{
      $ids_or_names = $self->set_names; 
      $method = 'fetch_by_name';
    }
  
    foreach my $var(@{$ids_or_names}){
      my $set = $set_adaptor->$method($var);
      
      if(! defined $set ){
        push @failed, $var;
      }
      elsif($handle_reps){
        
        if($set_type eq 'ResultSet'){
          #Validate there is just 1 supporting signal InputSubset (which is a replicate)
          
          #Check support type
          #Don't bother handling input_set here as it is being dropped
          if($set->table_name ne 'input_subset'){
            push @failed, $var;
            #$self->throw_no_retry('Expected input_subset support but found '.
            #  $set->table_name." support for ResultSet:\t".$set->name);  
          }
          else{
            
            my @issets = grep { ! $_->is_control } @{$set->get_support};
                   
            if(scalar(@issets) != 1 ||
               (! $issets[0]->replicate) ){
              #Replicate paranoia here
              #InputSubsets should be replicates
              #but historically we have had merged InputSubsets
              push @failed, $var;
            }    
          }   
        }
        elsif(! $set->replicate){ #Should be InputSet
          push @failed, $var;
        }
      }
      else{
        push @$sets, $set;
      }
    }
  }
  else{ #Must be other filters  
    my $constraints = $self->constraints_hash;
     
    #Need to account for analysis or format
    #i.e. we don't want to queue up the Segmentation input_sets
    #This should not really require and input_set
    #and should be loaded like and external set i.e. feature_set only 
   
    #Add string_param_exists here to validate states
    $sets = $set_adaptor->fetch_all( {constraints         => $constraints,
                                      string_param_exists => 1} ); 
    
    #Alternatively we could add a contraint for is_replicate
    #not for ResultSets at present, as they do not have the replicate field
                                        
    if($handle_reps){  
      my @rep_sets;
      
      foreach my $set(@$sets){
        
        if($set_type eq 'InputSet' &&
           $set->replicate){
          push @rep_sets, $set;  
        } 
        elsif($set_type eq 'ResultSet' && 
              ($set->table_name eq 'input_subset')){
          
          my @issets = grep { ! $_->is_control } @{$set->get_support};
                   
          if(scalar(@issets) == 1 && $issets[0]->replicate){
             push @rep_sets, $set;  
          }
        }
      }
      
      @$sets = @rep_sets;
    }
  }
    
       
  if(@failed){
    $throw = 'Failed to identify some '.$self->set_type." sets. Names or IDs don't exist";
    
    if($handle_reps){
      $throw .= ', and/or \'replicates\' were requested but could not find unique signal replicate'; 
    }
    $throw .= ":\n\t".join("n\t", @failed);
  } 
  
  
  ### FETCH/CACHE INPUTSUBSET CONTROLS ###
  @failed = ();
  my (@failed_ctrls, %ctrl_cache);
  my $x_grp_ctrls = $self->param_silent('allow_inter_group_controls');
  my $use_exp_id;
  
  if($set_type eq 'InputSubset'){
    
    if($self->control_experiments){
      my $exp_adaptor = $self->out_db->get_ExperimentAdaptor;
      my ($exp, $exp_group, $exp_method);
      
      #test whether we have a dbID or a name, assume first is representative of all
      if($self->control_experiments->[0] =~ /^[0-9]+$/){
        $exp_method = 'fetch_by_dbID';
      }
      else{#Assume we have a name
        $exp_method = 'fetch_by_name';
      }    
     
      foreach my $var(@{$self->control_experiments}){
        $exp = $exp_adaptor->$exp_method($var);
        
        if(! defined $exp){
          push @failed, $var;
          next; #$var
        }
        #This should really be Study! Not experimental group. So is not 100% safe
        $exp_group = $exp->get_ExperimentalGroup->name;
        
        $self->_cache_controls(
          $set_adaptor->fetch_all({constraints => {experiments => [$exp]}}),
          \%ctrl_cache, 
          \@failed_ctrls,
          undef,
          $x_grp_ctrls);
        #populates %ctrl_cache and @failed_ctrls
      }
     
      if(@failed){
        $throw .= "\nFailed to fetch some control_experiments. Names or IDs don't exist:\n\t".
                   join("n\t", @failed);   
      } 
      
      if(@failed_ctrls){
        $throw .= join("\n", @failed_ctrls); 
      }    
    }
     
    $use_exp_id = 1 if ! $self->identify_controls;
    #The subtlety here is that we may use a control for an isset which
    #which is not in the control_experiment and not directly associated
    #with that experiment. i.e. we are turning on identify_controls
    #for those experiments which aren't supported by control_experiments  
    #This is why we need to exp_name in the clabel!
    #We can simply test both keys below
    #Which will allow us to use control_experiments or experiment associated controls
    #If they are both defined, then they must match, else we don't know which ones to use!
    #Using associated controls would simply require omitting control_experiments
    #Over-riding associated controls with control_experiments would require deleting them 
    #from the DB
    #Or we could have a control over-ride flag?
    #Which would either use control_experiments preferentially or exclusively?
    
    #The other option is to simply turn on identify_controls, when we specify 
    #control_experiments?
    #This would cause problems when we have the simple case of having controls associated
    #with each experiment
        
    @failed_ctrls = ();
    $self->_cache_controls($sets,
                           \%ctrl_cache,   
                           \@failed_ctrls,
                           $use_exp_id,
                           $x_grp_ctrls);
    if(@failed_ctrls){
      $throw .= join("\n", @failed_ctrls);  
    }
  }

  
  ### SET UP VALIDATION VARS ###
  my ($no_rel_date, $rel_date, $force, $ignore, $rel_month, $tracking_adaptor);

  if( $self->validate_InputSubset_tracking ) { #we know this is an input_set
    $tracking_adaptor = $self->tracking_adaptor;
  
    $rel_month   = $self->get_param_method('release_month',    'silent');
    $no_rel_date = (defined $rel_month) ? 0 : 1;
    #currently in american format, hence we have a release month for safety
    #rather than a full date string  
    #Don't handle days here as release day is likely to shift, 
    #and we can just manage by reseeding with force_embargoed?
  }
  

  ### BUILD THE OUTPUT ID FOR EACH SET ###
  my (%output_ids, %seen_issets, %embargoed_issets, %to_download_issets, %control_reqd);
  my $peak_module;
  my $x_group_ctrls    = $self->param_silent('allow_inter_group_controls');
  my $allow_no_ctrls   = $self->allow_no_controls;
  my $force_embargoed  = $self->force_embargoed;
  my $ignore_embargoed = $self->ignore_embargoed;
  
  
 SET: foreach my $set( @$sets ){
    my $embargoed   = 0;
    my $to_download = 0;
  
      
  
  
    if( $self->validate_InputSubset_tracking ){ 
      my @vsets = ($set);
     
      if($set_type eq 'InputSet'){
        @vsets = $set->get_InputSubsets;  
      }
     
      foreach my $vset(@vsets){
       
        if(! exists $seen_issets{$vset->dbID}){
          $seen_issets{$vset->dbID} = undef;
  
          if( my @embargoed = $tracking_adaptor->is_InputSubset_embargoed($set, $rel_month) ){ 
   
            if(! ($force_embargoed || $ignore_embargoed)){
              my $rd_txt = '';
     
              if($no_rel_date){
                $rd_txt = "\nOr maybe you want to specify a release_month when seeding this analysis?"; 
              }
     
              $throw .= "\nFound InputSubset(s) which is not out of embargo:\n\t".
               join("\n\t", @embargoed).
               "\nYou can over-ride this by specifying force_embargo or ignore_embargo.".
               $rd_txt;
            }
            elsif($ignore_embargoed){
                         
              if($set_type eq 'InputSubset'){
                  $embargoed = 1;
                }
              }
              else{#input_set
                next SET;
              }
            }
          }
      
          if(! $tracking_adaptor->is_InputSubset_downloaded($set)){
            $to_download = 1;
            $throw .= "\nInputSet has InputSubsets which are not downloaded:\t\n".$set->name; 
            #'Please run download_input_set_data or specify allow_downloads
            #todo throw for now until we have written an analysis module
            #build the download output_id or mark it for batch flow to the download analysis
            #we shouldn't dataflow to DefineInputSet directly if we
            #are missing some fastqs!
            #Needs some more work, as if a control is missing, then we will need to batch up all the 
            #dependant InputSubset groups        
         }
         
         #Cache embargoed/to_download ctrls
     
         
         if($set->is_control &&
          ($embargoed || $to_download)){ #Set embargoed/to_download in the ctrl cache
           #Recreate the cache key logic here
           my $clabel   = $set->cell_type->name;
    
           if(! $x_group_ctrls){
             $clabel .= '_'.$set->get_Experiment->experimental_group_id;# <- probably doesn't exist (write wrapper and method!)
           }
        
           if(exists $ctrl_cache{$clabel} &&
              exists $ctrl_cache{$clabel}->[0]->{$set->dbID}){
              #double exists so we don't auto vivify 
              
              if($embargoed){
                $ctrl_cache{$clabel}->[3] ||= [];
                push @{$ctrl_cache{$clabel}->[3]}, $set->dbID;
              }
              
              if($to_download){
                $ctrl_cache{$clabel}->[4] ||= [];
                push @{$ctrl_cache{$clabel}->[4]}, $set->dbID;
              }
           }
            
                
           if($use_exp_id){
             $clabel   = $set->cell_type->name.'_'.$set->experiment_id;
           
             if(exists $ctrl_cache{$clabel} &&
                exists $ctrl_cache{$clabel}->[0]->{$set->dbID}){
                #double exists so we don't auto vivify 
              
               if($embargoed){
                 $ctrl_cache{$clabel}->[3] ||= [];
                 push @{$ctrl_cache{$clabel}->[3]}, $set->dbID;
               }
              
               if($to_download){
                 $ctrl_cache{$clabel}->[4] ||= [];
                 push @{$ctrl_cache{$clabel}->[4]}, $set->dbID;
               }
             }
           }
         }
       }
     }#end of validate_InputSubset_tracking
   
      
     if($set_type eq 'InputSubset'){
      
       if(! $set->is_control){    
         #We already have the ctrls in the ctrl_cache
        
         #automatically use control attached to the same experiment
         #how does this interact with control_experiments?
         #This should not cache as we don't want to identify_controls between experiments
         #too much? just auto identify controls?
        
         my $clabel   = $set->cell_type->name;
      
         if(! $x_group_ctrls){
           $clabel .= '_'.$set->get_Experiment->experimental_group_id;# <- probably doesn't exist (write wrapper and method!)
         }
        
         my $non_exp_ctrl = $ctrl_cache{$clabel} if exists $ctrl_cache{$clabel};
         my $exp_ctrl;
        
         if($use_exp_id){
           $clabel   = $set->cell_type->name.'_'.$set->experiment_id;
           $exp_ctrl = $ctrl_cache{$clabel} if exists $ctrl_cache{$clabel};
         }
        
         if( (! ($non_exp_ctrl || $exp_ctrl)) &&
            ! $allow_no_ctrls){
           $throw .= "Could not find a control for InputSubset:\t"     
         }
         elsif($non_exp_ctrl && $exp_ctrl){
           #Make sure these are the same!
          
           if($non_exp_ctrl->[2] ne $exp_ctrl->[2]){
             $throw .= "Found control clash between experiments:\t$non_exp_ctrl->[2] & $exp_ctrl->[2]";
           }
        
           #Can we assume these are identical?
           #They might not be control_experiments will get all the control from the experiments
           #but controls from the sets we have fetched (either identify_controls or just associated controls)
           #might not represent all the controls from an experiment
           #In this case, we can take the control_experiments
           #As this won't break any inter/intra experimental_group rules
           #Just do this below
         }
         else{
        
           #Preferentially take possibly more complete control_experiment derived controls
           $non_exp_ctrl ||= $exp_ctrl; 
           #Now we need to group these controls with their signal sets
           #as we need to flow these in a batch  
         
           #We need to decide on InputSet name here
           #experiment name or experiment_name with TR number
           #There is a change here that we have identified a subset of subsets for a given
           #experiment, so this may be slightly misleading
           #No easy way of resolving the naming for now, so let's just
           #leave it for now
    
           my $ctrl_cache_key;
        
           if(! $non_exp_ctrl){
             $ctrl_cache_key = 'no_controls';
             
             ### VALIDATE feature_set_analysis doesn't needs control (regardless of allow_no_controls )###
             #Here rather than in DefineResultSets to catch it earlier and prevent 
             #doomed job/input_id creation which may clash with future jobs submission and 
             #hence require some sort of recover/rollback flag
             
             #Do we need to be mindful of result_set_only here?
             #still don't want to create erroneous ResultSet even if we aren't running
             #the peak analysis. But wiggle track really doesn't require control
             #This is a small corner case, just bail out here for now.
             #But if we wanted to allow it would also need to test in DefineMergedOutputSet
             #Would also need to support rollback ResultSet to import with updated controls
             #this will result in collection being regenreated unecessarily :/
             
             #no_result_set is now redundant!!! As result_sets will replace input_sets
             #no_result_set behaviour will just change such that we don't create and associate a collection
             #and we don't set the relevant states
    
             ### GET FEATURE SET ANALYSIS ###
             my $fset_anal = $self->feature_set_analysis;
    
             if(! defined $fset_anal){
             
               if(exists $self->default_feature_set_analyses->{$set->feature_type->name}){
                 $fset_anal = $self->default_feature_set_analyses->{$set->feature_type->name}; 
               }
               elsif(exists $self->default_feature_set_analyses->{$set->feature_type->class}){
                 $fset_anal = $self->default_feature_set_analyses->{$set->feature_type->class};
               }
               else{
                 $throw .= "No default feature_set analysis available for ".$set->name.
                   ".\n Please add FeatureType name or class to  default_feature_set_analyses ".
                   "in BaseSequencingAnalysis::default_options or specify -feature_set_analysis\n"; 
                 next SET;
               }
             }
              
             ### TEST WHETHER IT NEEDS CONTROLS ###          
             if(! exists $control_reqd{$fset_anal}){
               eval { $fset_anal = &scalars_to_objects($self->out_db,
                                                       'Analysis', 
                                                       'fetch_by_logic_name',                                                        [$fset_anal])->[0]; };
               if( ! $fset_anal){
                 $throw .= 'Could not get analysis for '.$set->name."\n$@\n";  
                 next SET;
               }
            
               eval { $peak_module = $self->validate_package_from_path($fset_anal->module); };
              
               if(! $peak_module){
                 $throw .= 'Failed to validate analysis module for '.$set->name."\n$@\n"; 
                 next SET;                  
               }
                
               $control_reqd{$fset_anal} = $peak_module->requires_control;
             }
            
             if($control_reqd{$fset_anal}){
               $throw .= $set->name.' requires controls for '.$fset_anal->logic_name.
                 ", but none have been identified\n";
               next SET;
             }
           }
           else{        
             $ctrl_cache_key = $non_exp_ctrl->[0]->experiment->dbID;
           }
       
           
           ### BUILD OUTPUT ID ###
           
           if(! exists $output_ids{$ctrl_cache_key}){
             
             #Handle ignore_embargoed controls, and therefore all depedant subsets 
             
             $output_ids{$ctrl_cache_key} = {%$batch_params,
                                             input_subset_ids     => {},  #keys are set names or 'controls'
                                             %$dataflow_params}; 
             #Add ctrls first
             if($ctrl_cache_key ne 'no_controls'){
               #This should get updated with embargoed/to_download status
               #of subsequent ctrl sets which haven't yet been seen
               $output_ids{$ctrl_cache_key}->{input_subset_ids}->{controls} = $non_exp_ctrl;  
             }       
           }
  
           #Now add sig reps based on the exp/cell_type/feature_type input_subset_ids
           #Currently exp ID is a synecdoche for this, but let's future proof
           #in case we do move to proper experiment definition i.e. a study with
           #>1 ctype/ftypes           
           #my $iset_cache_key = $set->experiment->dbID.'_'.
           #                       $set->cell_type->dbID.'_'.
           #                       $set->feature_type->dbID;         
           #Cache key is essentially the set name
           #although this is more expensive to generate than a key just based on IDs
           
           #Define set name
           
           
           #my $ctype    = $set->cell_type->name;
           #my $ftype    = $set->feature_type->name;
           #my $set_name = $self->parse_study_name($set->experiment->name, 
           #                                        $set->cell_type->name,
           #                                        $set->feature_type->name);       
           #This is slightly odd as we are trying to future proof against
           #the experiment names changing to more generic study names
           #We add the ctype and ftype to form the iset name so we need
           #to strip these off if they are already present
             
           #Probably need to put this in a Base method           
           #$set_name = $ctype.'_'.$ftype.'_'.$set_name;
           
           #set name here will get appended with analysis logic name, but we ideally want 
           #the TR1 to be at the end, so need to append the logic_name first
           #we already have the alignment analysis, but this is really polluting
           #this module
           #we should just sub key on the rep number, and let DefineResultSets handle it.
           
           #This prefix should always be the standard (i.e. not the control) prefix
           my $set_name = $self->get_set_prefix_from_Set($set);        
           #my $rep      = $set->replicate;
           
           if($embargoed){
             $embargoed_issets{$ctrl_cache_key}{$set_name} ||= [];
             push @{$embargoed_issets{$ctrl_cache_key}{$set_name}}, $set->dbID;
           }
           
           if($to_download){
             #do we need to add dbids and replicate keys in here?
             $to_download_issets{$ctrl_cache_key}{$set_name} ||= [];
             push @{$to_download_issets{$ctrl_cache_key}{$set_name}}, $set->dbID;  
           }
           
           
           #why do we need to key on rep here?
           #can't we just have an arrayref value here?
           
           $output_ids{$ctrl_cache_key}->{input_subset_ids}->{$set_name} ||= [];
           push @{$output_ids{$ctrl_cache_key}->{input_subset_ids}->{$set_name}}, $set->dbID;
         }        
       }
     }
     else{ #Not an input_subset (currently just a ResultSet)
      my $group_key = 'ungrouped';
     
      if($self->only_replicates){
        #test is replicate here?
        #probably already tested previously     
        
        #just parent ResultSet name (which doesn't exist yet)
        #This is a little fragile as it depends on our nomenclature
        #would replace this with the get set prefix
        #then at least we are using the same code everywhere
        #($group_key = $set->name) =~ s/_TR[0-9]+$//o;
        
        $group_key = $self->get_set_prefix_from_Set($set).'_'.$set->analysis->log_name;    
        #Appened analysis name here, just in case we are running
        #multiple analyses for the same set of data.
      }
     
    
       $output_ids{$group_key} ||= [];
       push @{$output_ids{$group_key}}, {%$batch_params,
                                         %$dataflow_params,
                                         dbID     => $set->dbID, 
                                         set_name => $set->name};
    }
 

    #Need to submit input_subsets them based on shared controls!
    #As we need to wait for the control jobs to download/align 
    #first before we can perform any other analyses?
    

    #Is this true?
    #This is only true if the control has not already been aligned!
    #as we need to semaphore downstream analysis based on success of control alignment
    #this does not however stop us from doing the signal alignments
    #we just don't want to sempahore the whole fan, as we don't want jobs blocking each other
    #only the control job
    
    #WE can't actually do this as this would require semaphore from different analysis
    
    #Can we do this with an accu output from the control?
    #No, as the control job may not finish before the signal jobs dataflow
    #Also, this would create failures of downstream analyses, when really on the 
    #control analysis has failed.
    #This is messy and also pollutes the error output, so real errors would be easily overlooked
  } 
    
  
  #Now iterate through output_ids doing the right thing
  
  my $warn_msg = '';
  
  #This needs completely changing as the keys and values have now changed!!
  
  
 KEY: foreach my $key(keys %output_ids){
    #we need the key to access embargoes/to_download_issets cache 
    my $oid = $output_ids{$key};  
    
    #Also need to check the embargoed status, and skip InputSets/Subsets which are dependant on 
    #an embargoed subset
    #no write to support listing sets before actually seeding them
    #need to be able to run this as Stand alone job, but with access
    #to config. Leo is on the case here. 
    
    if($set_type eq 'InputSubset'){
      #{%$batch_params,
      # input_subset_ids     => {},  #keys are set names or 'controls'
      # %$dataflow_params}; 
 
      my $branch = 3; #ignore or not embargoed 
      
      if(exists $oid->{input_subset_ids}->{controls}){
        
         #[{$cexp_isset->dbID => $cexp_isset, ...},
         #                       $cexp_isset->feature_type->name,
         #                       $exp_id,
         #                       [],  #embargoed sets
         #                       []   #to_download sets
         #        ];  
        
        
        if(scalar(@{$oid->{input_subset_ids}->{controls}->[3]}) > 0){ #CTRLs are embargoed!
          #build input_set specific log msg
          #appending to appropriate throw/warn string and setting branch
          my $ctrls = delete $oid->{input_subset_ids}->{controls};
          my $msg   = "\nFound embargoed control subsets ($key):\t".
                        join(', ', @{$ctrls->[3]})."\nRelated InputSets:\n\t".
                        join("\n\t", keys %{$oid->{input_subset_ids}})."\n";
          
          if($ignore_embargoed){
            $warn_msg .= "Skipping dataflow for the following InputSets.\n$msg";   
            $branch = 0;            
            next KEY; #Not interested in InputSet specific embargoes
          }
          elsif($force_embargoed){
            $warn_msg .= "Forcing dataflow for embargoed controls ($key):\t".
                          join(', ', @{$ctrls->[3]})."\n";
            $oid->{input_subset_ids}->{controls} = $ctrls;   #add ctrls back in!
          }
          else{
            $throw .= $msg;
            $branch = 0;
            #don't next as we are probably interested in InputSet specific embargoes?
          }
        }         
      }
    
          
     #todo The more helpful messages about release_month and force/ignore_embargoed are above
     #we should print them once at the bottom
      
      if(exists $embargoed_issets{$key}){ #Some of the signal subsets are embargoed
        
        foreach my $iset_name(values %{$embargoed_issets{$key}}){
        
          if($ignore_embargoed){
            $warn_msg .= "\nSkipping dataflow for InputSet with embargoed signal subsets:\t$iset_name\n";   
            delete $oid->{input_subset_ids}->{$iset_name};
          
          }
          elsif($force_embargoed){
            $warn_msg .= "Forcing dataflow for InputSet with embargoed signal subsets:\t$iset_name\n";
          }
          else{
            $throw .= "Found InputSet with embargoed signal subsets:\t$iset_name\n";            
            $branch = 0;
          }
        }
      }
        
      #todo Test %to_download_issets, and fan these to $branch to 2 if $branch is true 
      #(i.e. it's not embargoed)
      #Or do we want to download regardless of embargo status?
      #Likely yes, but leave that for now, as figuring out config
      #Also flow to branch 2 aswell after we have fanned
      #Do this by building a job group and using cache_job_group
              
      if($branch){
    
        if($no_write){
          my $txt = 'InputSets identified (';
          
          if(exists $oid->{input_subset_ids}->{controls}){
            $txt .= 'shared controls';
          }
          else{
            $txt .= 'no controls';
          } 

          print STDOUT $txt."):\n\t".join("\n\t", (keys %{$oid->{input_subset_ids}}))."\n";
        }
        elsif(scalar(keys %{$oid->{input_subset_ids}}) > 1 ||
              ((scalar(keys %{$oid->{input_subset_ids}}) == 1) && 
               ! exists $oid->{input_subset_ids}->{controls} ) ){
          #$self->branch_output_id($oid, undef, $branch);   
       
         #Tidy up control cache
         #probably need to revise this when we add branching for download analyses
          
         if(exists $oid->{input_subset_ids}->{controls}){
           $oid->{input_subset_ids}->{controls} = [keys(%{$oid->{input_subset_ids}->{controls}->[0]})];      
         }
         
            
          $self->branch_job_group($branch, [$oid]);              
        }
        else{#We have no signal subsets to flow!
          $warn_msg .= "No input_subsets to dataflow for control group:\t$key\n";
        }
      }
    }
    else{ #Not input_subset
     
      if($no_write){
        print STDOUT "\t".$oid->{set_name}.' ( '.$oid->{dbID}." )\n";
      }
      else{
     
        #Branches now generic wrt set type   
        #2 group of sets (optionally grouped by either controls or parent sets semaphored by 3)
        #>3 fans of those grouped sets, which require different processing
        #actually output id would change dependant on the set type, and the grouping strategy
        #grouping would have to be defined by analysis parameter
        #This will allow expansion of branches without renumbering
        #As branch 2 is always the semaphore, we should always flow this last
        
        #This does however mean that for dataflow with no semaphore we will never define
        #a funnel, and therefore only flow down branch 3 and not branch 2.
        
        my @job_group = (3, $oid);
        
        if($key ne 'ungrouped'){
          #Now define funnel job id with all set details i.e. RunIDR
          
          my $funnel_oid = {%$batch_params,
                            %$dataflow_params,
                            dbIDs     => [], 
                            set_names => []}; 
                            
          foreach my $id(@$oid){
            push @{$funnel_oid->{dbIDs}},     $id->{dbID};
            push @{$funnel_oid->{set_names}}, $id->{set_name};
          }
          
          push @job_group, (2, [$funnel_oid]);
        }
        
        
        $self->branch_job_group(@job_group);  
        
        #I don't think we are using this branch_config anymore??
        #if(defined $self->branch_config){
        #  #branch_key_method will have been validated by init_branch_config   
        #  #Can't do this out of the loop unless we test for config defined
        #  #Is this used any more here?
        #  my $branch_key_method = $self->branch_key_method;
        #  $self->branch_output_id($oid, $self->$branch_key_method);  
        #}
      }
    }
  } #end foreach $oid
    
  #Only throw here rather than above so we get all the information
  #about what has failed  
  
  warn $warn_msg if $warn_msg;
  $self->throw_now_retry($throw) if $throw;  
  return;
}


sub _cache_controls{
  my ($self, $cexp_issets, $ctrl_cache, $failed, $use_exp_id, $x_grp_ctrls) = @_;    
  my $clabel; 
   
  foreach my $cexp_isset(@{$cexp_issets}){
       
    if($cexp_isset->is_control){   
      $clabel    = $cexp_isset->cell_type_name;
      my $exp_id = $cexp_isset->experiment_id;;
    
      if($use_exp_id){
        $clabel .= '_'.$exp_id;
      }
      elsif(! $x_grp_ctrls){
        $clabel .= '_'.$cexp_isset->get_Experiment->experimental_group_id;# <- probably doesn't exist (write wrapper and method!)
      }
    
      if(! exists $ctrl_cache->{$clabel}){
        #Hash required here to avoid redundant controls
        $ctrl_cache->{$clabel} = [{$cexp_isset->dbID => $cexp_isset},
                                  $cexp_isset->feature_type->name, #Used for validation
                                  $exp_id, #Used for validation
                                  [],  #embargoed sets
                                  []   #to_download sets
                                  ];  
      }
      elsif( ($cexp_isset->feature_type->name eq $ctrl_cache->{$clabel}->[1]) &&
             ($exp_id eq $ctrl_cache->{$clabel}->[2]) ){
        $ctrl_cache->{$clabel}->[0]->{$cexp_isset->dbID} = $cexp_isset;
      }
      else{
        #May need to improve this message
        push @$failed, $clabel.' control '.$cexp_isset->name." clashes with cache:\n\t".
          $ctrl_cache->{$clabel}->[1]."\n\t".$ctrl_cache->{$clabel}->[2];    
      }
    }#fi is_control
    
  }#end foreach $isset
  
  return;# $failed; 
}



#Todo 
#Enable a preview of what we are going to dataflow here
#This will be done using standaloneJob.pl once this can access config from the DB
#will probably have to add a no_data_flow flag, which will print out instead of
#flow the output_ids

sub write_output {  # Create the relevant jobs
  my $self = shift;  
  $self->dataflow_job_groups;
  return;
}


#Could move this to the TrackingAdaptor
#if the TrackingAdaptor also updated it's cache after download
#is not implemented above

sub is_InputSet_downloaded {
  my $self = shift;
  my $iset = shift;
  
  if($iset){
    #This method handles/validates InputSets too
    $self->{input_set_downloaded} = $self->tracking_adaptor->is_InputSubset_downloaded($iset); 
  }   
  
  return $self->{input_set_downloaded};
}

1;
