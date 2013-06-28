=pod 

=head1 NAME

Bio::EnsEMBL::Funcgen::Hive::IdenetifySetInputs;

=head1 DESCRIPTION

This module simply takes a list of paramters and a 'set_type' to identify 
Sets used as inputs for various parts of the analysis pipeline.

=cut

package Bio::EnsEMBL::Funcgen::Hive::IdentifySetInputs;

use base ('Bio::EnsEMBL::Funcgen::Hive::BaseDB');



#Bio::EnsEMBL::Funcgen::Hive(::Config)
#We don't need to discriminate between Runnables and RunnableDBs anymore
#Just name the modules accordingly!


use warnings;
use strict;
#use Bio::EnsEMBL::DBSQL::DBAdaptor;
#use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor; 
#use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw (strip_param_args generate_slices_from_names 
#                                               strip_param_flags run_system_cmd);
use Bio::EnsEMBL::Utils::Exception qw (throw);


#todo -slice_import_status?


my %set_adaptor_methods = 
(
  input_set   => 'get_InputSetAdaptor',
  result_set  => 'get_ResultSetAdaptor',
  feature_set => 'get_FeatureSetAdaptor',
  #data_set?
);


#Accessors to catch compile time errors and avoid uncaught typos
sub constraints_hash {    return $_[0]->param('constraints_hash');    }
sub set_type {            return $_[0]->param_required('set_type');   }
sub feature_types {       return $_[0]->param('feature_types');       }
sub cell_types {          return $_[0]->param('cell_types');          }
sub experimental_groups { return $_[0]->param('experimental_groups'); }
#sub input_sets {          return $_[0]->param('input_sets');         }
sub analyses {            return $_[0]->param('analyses');            }
sub set_ids {             return $_[0]->param('set_ids');             }
sub set_adaptor {         return $_[0]->param('set_adaptor');         }
sub set_names {           return $_[0]->param('set_names');           }
sub sets {                return $_[0]->param('sets');                }
sub states {              return $_[0]->param('states');              } 



sub fetch_input {   # fetch parameters...
  my $self = shift @_;
  $self->SUPER::fetch_input;
  
  
  #Validate set_type
  my $set_type = $self->set_type; 
  
  if(! exists $set_adaptor_methods{$set_type}){
    throw("The -set_type $set_type is not supported by IdentifySetInputs.\n".
          "Valid options are ".join(' ', keys(%set_adaptor_methods)) );      
  }
  else{
    #Set the method name value to reset the value as the actual adaptor
    my $method = $set_adaptor_methods{$set_type};
    $self->param('set_adaptor', $self->out_db->$method);
  }
  
  
  $self->process_params([qw(set_names set_ids)], 1);#optional
  
  #Parse comma separated lists of filters into arrayrefs of string or objects
  $self->param('constraints_hash',
               $self->process_params([qw(feature_types cell_types states
                                       analyses experimental_groups)],
                                       1, 1)  #optional/as array flags
               );      
  

  #Catch mutally exclusive filter params  
  
  if($self->feature_types ||
     $self->cell_types ||
     $self->experimental_groups ||
     $self->analyses ||
     $self->states){ 
     #all these are OR filters except states which is an AND filter
      
    if($self->set_names || $self->set_ids){
      throw('You have specified mutually exclusive filter params for the '.
            "IdentifySetInputs analysis\nPlease specify restrict to ".
            '-set_name or -set_ids or a combination other filters '.
            '(e.g. -experimental_groups -feature_types -cell_types -analyses -states');
    }
  }
  elsif(! ($self->set_names || $self->set_ids)){
    throw('You must specifiy some IdentifySetInputs fitler params either '.
          '-input_sets or -set_ids or a combination of '.
          '-feature_types -cell_types -experimental_groups -states -analyses');
  }
  elsif($self->set_names && $self->set_ids){
    throw('You have specified mutually exclusive filter params for the '.
            "IdentifySetInputs analysis\nPlease specify restrict to ".
            'set_names or -set_ids or a combination other filters '.
            '(e.g. -experimental_groups -feature_types -cell_types -analyses -states');  
  }


  return;
}



sub run {   # Check parameters and do appropriate database/file operations... 
  my $self = shift @_;
  my %sets;
  my $set_adaptor = $self->set_adaptor;
  my $throw = 0;
  
  #For set_ids and set_names, catch undef return types
  if($self->set_ids){
  
    foreach my $id(@{$self->set_ids}){
      my $set = $set_adaptor->fetch_by_dbID($id);
      $sets{$id} = $set;
      $throw = 1 if ! $set;
    }
  }
  elsif($self->set_names){

    foreach my $name(@{$self->set_names}){
      my $set = $set_adaptor->fetch_by_name($name);
      $sets{$name} = $set;
      $throw = 1 if ! $set;
    }    
  }
  else{ #Must be other filters  
    my $constraints = $self->constraints_hash;

    #$self->param('input_sets', $constraints);
    #$self->process_params($set_type.'s');
    #or add $self->set_and_process_params?


    if(($self->set_type eq 'data_set') || 
       ($self->set_type eq 'result_set'  )){   
      throw($self->set_type.
        ' adaptor does not yet support fetch_all({constraints => $constraints,})');
     }
     
    #Need to account for analysis or format
    #i.e. we don't want to queue up the Segmentation input_sets
    #This should not really require and input_set
    #and should be loaded like and external set i.e. feature_set only 
   
    #Add string_param_exists here to validate states
     
    foreach my $set( @{$set_adaptor->fetch_all( {constraints => $constraints,
                                                 string_param_exists => 1} )} ){ 
      #warn "Found set ".$set->name;     
      $sets{$set->dbID} = $set;
      
      #$self->helper->debug(1, $set->name.'( '.$set->dbID.")");
      
      if($self->param('no_write')){
          print STDOUT "\t".$set->name.'( '.$set->dbID.")\n";
      }
    }
  }
    
      
  if($throw){
    throw('Failed to fetch some '.$self->set_type." Sets using names or IDs:".
          join("n\t", (map {$_.' => '.$sets{$_}} keys %sets)));    
  }
  
  $self->param('sets', [values %sets]);
  return;
}



#Todo 
#Enable a preview of what we are going to dataflow here
#This will be done using standaloneJob.pl once this can access config from the DB
#will probably have to add a no_data_flow flag, which will print out instead of
#flow the output_ids

sub write_output {  # Create the relevant jobs
  my $self = $_[0];
  
  foreach my $set(@{$self->sets}){
    #Flow dbID and name for readability in hive DB
    #flows to batch jobs (branch 2) e.g. StoreRollbackSets    
    $self->dataflow_output_id({dbID => $set->dbID, set_name => $set->name}, 2);
          
    #This should return the job_id
    #if absent, then this is already present, and nothign would be returned
    #do we need to warn about this?
  }
  
  return;
}


1;
