=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Funcgen::Hive::DefineDataSet

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::Hive::DefineDataSet;

use warnings;
use strict;
 
use Bio::EnsEMBL::Utils::Exception qw (throw);
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(scalars_to_objects);

use base qw(Bio::EnsEMBL::Funcgen::Hive::BaseDB);

#This assumes the InputSet has been previously registered, 
#and now we want simply to define/fetch the data set, feature and result set based on these data.

#params
#name
#result_set_only
#These allow over-ride of defaults
#feature_set_analysis
#result_set_analysis
#default_feature_set_analyses (shouldn't alter)
#default_result_set_analyses (shouldn't alter)

#todo
#-slice_import_status?
#This now needs to take lists of dbIDs for InputSubsets
#These should all share the same control
#The InputSets will be created iteratively
#but the data flow will be done in a batch fashion, passing on sets of InputSet IDs
#as appropriate
#change this analysis to DefineSets? As it deals with InputSets, FeatureSets, ResultSet and DataSets
#or should we just create another analysis which is DefineInputSets?
#as it is almost entirely different?


sub fetch_input {   # fetch parameters...
  my $self = shift @_;
  $self->check_analysis_can_run; #Could put this in Base::fetch_input?
  $self->SUPER::fetch_input;  
  my $set = $self->fetch_Set_input('ResultSet');#defines result_set method

  #refactor this default_analysis method in BaseDB? As this is reused in other runnables
  my %default_analyses;
  my @set_types = ();# ('result_set');
  push @set_types, 'feature_set' if ! $self->get_param_method('result_set_only', 'silent');
 
  #TODO Validate default set analysis keys exist as feature_type class or name?
  #This would fail for species with low coverage i.e. some names may be absent


  #Can we move some of this into BaseSequenceAnalysis as it will need to be used
  #by Run_QC_and_ALigner to enable selective data flow into the relevant 
  #IdentifyInputSets analyses
  
  #we could simplify this loop now, as we will always have the result set
  #but don't we want to keep it generic, such that we can re-use it however we want?
  #This is already hardocded to take a result set there is currently no other possibility
  #apart from input_set, which will disappear
  
  
  #Do we even want to define a DataSet if there is no feature_set? 
  #durely we shoudl flow the relevant set types down each branch of the config
  #rather than just the data set?
  
  foreach my $set_type( @set_types ){  
    my $anal_type = $self->param_required($set_type.'_analysis_type'); 
    my $set_lname = $self->param_silent($anal_type.'_analysis');
    
    if(! defined $set_lname){
      
      $default_analyses{$set_type} = $self->param_silent('default_'.$set_type.'_analyses');     
      
      if(! defined $default_analyses{$set_type}){
        throw("Please define -${set_type}_analysis or add to default_${set_type}_analyses".
          " in the default_options config");
      }
      
     
      if(exists $default_analyses{$set_type}{$set->feature_type->name}){
        $set_lname = $default_analyses{$set_type}{$set->feature_type->name}; 
      }
      elsif(exists $default_analyses{$set_type}{$set->feature_type->class}){
        $set_lname = $default_analyses{$set_type}{$set->feature_type->class};
      }
      else{
        throw("No default $set_type analysis available for ".$set->feature_type->name.
          '('.$set->feature_type->class.
          ").\n Please add FeatureType name or class to  default_${set_type}_analyses ".
          "in default_options config or specify -${set_type}_analysis"); 
      }
    }                                  
    
    #Catch undefs in config
    if(! defined $set_lname){
      throw("Unable to identify defined $set_type analysis in config for ".$set->feature_type->name.
      ".\nPlease define in default_${set_type}_analyses config or specify -${set_type}_analysis");     
    }    
    #can't use process_params here as the param name does match the object name
    #We could over-ride this with a second arrayref of class names                                      
    $self->param($set_type.'_analysis', 
                 &scalars_to_objects($self->out_db,
                                     'Analysis', 
                                     'fetch_by_logic_name',
                                     [$set_lname])->[0]);                                              
  }
  

  return;
}


#TODO handle ResultSet only run
#This would take a flag and create the ResultSet in isolation
#would need to validate not FeatureSet stuff was set
#this might clash with some defaults if we had mixed set types in the pipeline
#i.e. some with and some without datasets
#feature_set_analysis is always created dynamically anyway, so this would be dependant on
#-result_set_only 


#TODO Add support for feature_set_only! IDR sets!


sub run {   # Check parameters and do appropriate database/file operations... 
  my $self   = shift;
  my $helper = $self->helper;  
  my $rset   = $self->ResultSet;
  my $set;    

  #todo migrate this to the Importer as define_OutputSet
  #there is some overlap here of param validation between BaseImporter and hive
  
      
  #Never set -FULL_DELETE here!
  #It is unwise to do this in a pipeline and should be handled
  #on a case by case basis using a separate rollback script    
  
  #Should also never really specify recover here either?
  #This bascailly ignores the fact that a ResultSet may be lin ked to other DataSets
  
  
  if( $self->result_set_only ){
    throw('Pipeline does not yet support creation of a ResultSet without and associated Feature/DataSet');
     
    #$set = $helper->define_ResultSet
    #  ( 
    #   -NAME                 => $rset->name,#.'_'.$rset_anal->logic_name,
    #   #-FEATURE_CLASS        => result | dna_methylation',
    #   #currently set dynamically in define_ResultSet
    #   -SUPPORTING_SETS      => [$rset],
    #   -DBADAPTOR            => $self->out_db,
    #   -RESULT_SET_ANALYSIS  => $self->param('result_set_analysis'),
    #   -RESULT_SET_MODE      => $self->param('result_set_mode'),
    #   -ROLLBACK             => $self->param('rollback'),
    #   -RECOVER              => $self->param('recover'),
    #   -SLICES               => $self->slices,
    #   -CELL_TYPE            => $iset->cell_type,
    #   -FEATURE_TYPE         => $iset->feature_type,
    # );
        
  }
  else{
    my $fset_anal = $self->param('feature_set_analysis');
    
    $set = $helper->define_DataSet
      (
       -NAME                 => $set->name.'_'.$fset_anal->logic_name,
       -FEATURE_CLASS        => 'annotated', #Is there overlap with rset feature_class here?
       -SUPPORTING_SETS      => [$rset],
       -DBADAPTOR            => $self->out_db,
       
       -FEATURE_SET_ANALYSIS => $fset_anal,
       -RESULT_SET_ANALYSIS  => $self->param('result_set_analysis'),
       -RESULT_SET_MODE      => $self->param('result_set_mode'),
       -ROLLBACK             => $self->param('rollback'),
       -RECOVER              => $self->param('recover'),
       -SLICES               => $self->slices,    
       -CELL_TYPE            => $rset->cell_type,
       -FEATURE_TYPE         => $rset->feature_type,
       #-DESCRIPTION   => ?
       #-DISPLAY_LABEL => ?
      );  
  }


  #No tracking required here?
  #Todo review whether rollback handles status entries correctly
  
  
  
  #TODO Need to implement branch config here!
  

  #Add set_type here as result_set_only could be change between writing
  #this output_id and running a down stream analysis  
  $self->param('output_id', 
               {%{$self->batch_params},
                dbID       => $set->dbID, 
                set_name   => $set->name,
                set_type   => ($self->result_set_only ? 
                               'ResultSet' : 'DataSet')}
              );
  
  return 1;
}


sub write_output {  # Create the relevant jobs
  my $self = shift;
  $self->helper->debug(1, 'DefineDataSet data flowing:', $self->param('output_id'));
  $self->dataflow_output_id($self->param('output_id'), 1);
  return;
}


1;
