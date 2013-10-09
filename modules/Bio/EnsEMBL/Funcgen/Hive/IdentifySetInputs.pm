=pod 

=head1 NAME

Bio::EnsEMBL::Funcgen::Hive::IdenetifySetInputs;

=head1 DESCRIPTION

This module simply takes a list of paramters and a 'set_type' to identify 
Sets used as inputs for various parts of the analysis pipeline.

=cut

package Bio::EnsEMBL::Funcgen::Hive::IdentifySetInputs;

use warnings;
use strict;

use Bio::EnsEMBL::Utils::Exception qw (throw);
use base ('Bio::EnsEMBL::Funcgen::Hive::BaseDB');



#todo -slice_import_status?


my %set_adaptor_methods = 
(
  input_subset => 'get_InputSubsetAdaptor',
  input_set    => 'get_InputSetAdaptor',
  result_set   => 'get_ResultSetAdaptor',
  feature_set  => 'get_FeatureSetAdaptor',
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



sub fetch_input {   # fetch parameters...
  my $self = $_[0];
  $self->SUPER::fetch_input;
  
  #$self->helper->debug(3, "IdentifySetInput params:\t", $self->input_job->{_param_hash});
  
  
  #Validate set_type
  my $set_type = $self->get_param_method('set_type', 'required'); 
  
  if( $self->get_param_method('validate_InputSubset_tracking', 'silent') ){
    
    if($set_type ne 'input_set') {
      throw('It is not possible to validate_embargo for '.$set_type.
        's. The IdentifySetInputs analysis is misconfigured. Please omit validate_embargo '.
        'or correct the set_type to \'input_set\'');
    }
    
    if($self->set_param_method('force_embargoed',  'silent') &&
       $self->set_param_method('ignore_embargoed', 'silent') ){
      throw('force_embargoed and ignore_embargoed are mutually exclusive parameters. '.
        'Please omit one or both');     
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
  
  
  $self->process_params([qw(set_names set_ids replicates)], 1);#optional
  
  #Parse comma separated lists of filters into arrayrefs of string or objects
  $self->set_param_method('constraints_hash',
                          $self->process_params([qw(feature_types cell_types states
                                                    analyses experiments experimental_groups)],
                                                1, 1)  #optional/as array flags
                          );      
  

  #Catch mutally exclusive filter params  
  
  #Add experiment support here for InputSubsets
  #InputSubset adaptor only supports experiments, names and ids
  #Or can we just do this by query composition to the experiment table?
  
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
            '(e.g. -experimental_groups -exoeriments -feature_types -cell_types -analyses -states');
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


  $self->init_branch_config(1, 1);
  #Flags are:
  #Optional branch config, as we only need this for IdentifyAlignInputsets 
  #Validate branch_key_method if we do have config

  return;
}



sub run {   # Check parameters and do appropriate database/file operations... 
  my $self = shift @_;
  my $sets = [];
  my $set_adaptor = $self->set_adaptor;
  my $throw = 0;

  #Can grab these here as know these params will not change as we iterate
  #through sets below
  #These will only be flowed to the next job
  my $dataflow_params = $self->dataflow_params(1);#optional flag
  #where as batch_params flow across all jobs in this seeded batch
  #(for those that support/require batch params)
  my $batch_params    = $self->batch_params; 
  my (@failed_sets);
  my $reps = $self->replicates;
  
  #Can't batch flow replicates, as this would cause 
  #the following to fail when dealing with other non-replicate sets
  #replicates should be an analysis_param  
  
  #For set_ids and set_names, catch undef return types
  if($self->set_ids){
  
    foreach my $id(@{$self->set_ids}){
      my $set = $set_adaptor->fetch_by_dbID($id);
      
      if(! defined $set ||
         ($reps && (! $iset->replicate)){
        
        push @failed_sets, $id;
        $throw = 1;
      }else{
        push @$sets, $set;
      }
    }
  }
  elsif($self->set_names){

    foreach my $name(@{$self->set_names}){
      my $set = $set_adaptor->fetch_by_name($name);
        
      if(! defined $set ||
         ($reps && (! $iset->replicate)){
        
        push @failed_sets, $name;
        $throw = 1;
      }else{
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
    
    #We could add a contraint for is_replicate
                                      
    if($reps){  
      foreach my $set(@$sets){
        
        if(! $set->replicate){
          push @failed_sets, $set->name;  
        } 
      }
    }
  }
    
       
  if($throw){
    $throw = 'Failed to identify some '.$self->set_type." sets. Names or IDs don't exist";
    
    if($reps){
      $throw .= ', and/or they are merged sets and \'replicates\' were requested'; 
    }
    
    throw($throw.':'.join("n\t", @failed_sets));    
  } 
  
  
  #Here we need to do the InputSubset embargo date stuff
  #but this is only for IdentifyAlignInputSets
  #so we need a validate_embargo flag as default for this analysis
  #removing this would essentially be the same as setting force_embargoed
  #but we should hardcode validate_embargo in the config
  #and force the usage of ignore|force_embargoed
  
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
  
  if($self->param_silent('no_write')){
    print STDOUT "Identified the following ".$self->set_type."s:\n";  
  }

#todo Could do with a way of listing embargoed
#when no_write is defined
  
 SET: foreach my $set( @$sets ){
         
    #This need to be a more general flag
    #validate_InputSubset_tracking
         
         
    if( $self->validate_InputSubset_tracking ){ #we know this is an input_set
    
      if( my @embargoed = $tracking_adaptor->is_InputSubset_embargoed($set, $rel_month) ){ 
     
        if(! ($self->force_embargoed || $self->ignore_embargoed)){
       
          my $rd_txt = '';
       
          if($no_rel_date){
            $rd_txt = "\nOr maybe you want to specify a release_month when seeding this analysis?"; 
          }
       
          throw("Found InputSubset(s) which is not out of embargo:\n\t".
            join("\n\t", @embargoed).
            "\nYou can over-ride this by specifying force_embargo or ignore_embargo.".
            $rd_txt);
        }
        elsif($self->ignore_embargoed){
          next SET; 
        }
      }
        
        
      #Dang, we already have branch_config for this
      #Need to hardcode this as the branch
       
    
      #if($branch_key_method eq 'InputSubset_downloaded'){
    
          #throw for now until we have written an analysis module
          #to handle the downloading      
      if(! $tracking_adaptor->is_InputSubset_downloaded($set)){
        throw("InputSet has InputSubsets which are not downloaded:\t\n".$set->name); 
         
          #'Please run download_input_set_data or specify allow_downloads
              
                   
      }
    }
    
    
    my $output_id = {%$batch_params,
                     dbID => $set->dbID, 
                     set_name => $set->name,
                     %$dataflow_params};
       
    #no write to support listing sets before actually seeding them
    #need to be able to run this as Stand alone job, but with access
    #to config. Leo is on the case here. 
    if($self->param_silent('no_write')){
      print STDOUT "\t".$set->name.' ( '.$set->dbID." )\n";
    }
    else{
      $self->branch_output_id($output_id, undef, 2);
        
      #handle other branch here too (generically)
      #Will only ever have 1 per instance, but this is depedant on 
      #alignment_analysis for IdentifyAlignInputSets
      if(defined $self->branch_config){
        
        #branch_key_method will have been validated by init_branch_config   
        
        #really need to be able to pass some params to this method
        #but no way of doing this generically as the method could be anything!
        
        #we want to be able branch based on the download state of an InputSubset
        #basically flow to an intermediate download analysis 
        #before we flow to the align analysis
        
        #Can't do this out of the loop unless we test for config defined
        my $branch_key_method = $self->branch_key_method;
           
        $self->branch_output_id($output_id, $self->$branch_key_method);  
      }
    }
  }  
    
  return;
}



#Todo 
#Enable a preview of what we are going to dataflow here
#This will be done using standaloneJob.pl once this can access config from the DB
#will probably have to add a no_data_flow flag, which will print out instead of
#flow the output_ids

sub write_output {  # Create the relevant jobs
  my $self = $_[0];  
  $self->dataflow_branch_output_ids;
  return;
}

sub is_InputSet_downloaded {
  my ($self, $iset) = @_;
  
  if($iset){
    #This method handles/validates InputSets too
    $self->{input_set_downloaded} = $self->tracking_adaptor->is_InputSubset_downloaded($iset); 
  }   
  
  return $self->{input_set_downloaded};
}

1;
