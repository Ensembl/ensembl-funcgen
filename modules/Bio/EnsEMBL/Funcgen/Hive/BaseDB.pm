=pod 

=head1 NAME

Bio::EnsEMBL::Funcgen::RunnableDB::Funcgen

=head1 DESCRIPTION

'Funcgen' is a base class for other runnables of the Funcgen Hive Pipeline
It performs common tasks such as connecting to the EFG DB etc...

=cut

package Bio::EnsEMBL::Funcgen::Hive::BaseDB;

use warnings;
use strict;

use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor; 
use Bio::EnsEMBL::Funcgen::DBSQL::TrackingAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( generate_slices_from_names assert_refs );
use Bio::EnsEMBL::Utils::Exception         qw( throw warning );

use base ('Bio::EnsEMBL::Funcgen::Hive::Base');


#This defines a set of parameters based on given parameters from the pipeline:
#All params need to be in pipeline_wide_params else we will get warnings e.g.
#ParamWarning: value for param('output_dir') is used before having been initialized!
#However we don't really want these to be pipelinewide!
#i.e. we don't want 

#Is the solution here to remove all the defaults where possible, if we are not passing them
#as pipeline wide or default analysis param
#seed_pipeline or dataflow with add these
#but there is still no way around the ParamWarnings

#Need to be mindful about flowing things like slices etc.


#params
#db_output_dir
#slices
#skip_slices

sub fetch_input {
  my $self = $_[0];
  $self->SUPER::fetch_input;


  
  #Set the DB if we haven't already in Base (for tracking)
  #This also sets species and assembly methods
  $self->_set_out_db if ! $self->use_tracking_db;
    
 
  #Do not force mandatory dir params heres
  #let subclasses do this using validate_dir_param
  #this is because not all otuptu will go in the same dir
  #i.e. alignment may use BaseDB or will it?
  #no it will use Base with a tracking DB
  #BaseDB is really where the output it the DB
  
  #Set the base output_dir for analyses which output to the DB
  
  $self->validate_dir_param
   ('db_output_dir',
    undef,  #create/writeable? #should this be undef, and force the runnables to validate_dir_param output_dir?
    $self->data_root_dir.'/output/'.$self->out_db->dbc->dbname); #default
  
  
  $self->process_params(['slices', 'skip_slices'], 1, 1);#optional/as array flags
 
  #validate/generate slices/skip_slices here
  if($self->slices || $self->skip_slices){
    $self->slice_objects; 
    #Validates custom slice & skip slices before we start fanning jobs
  }

  return;
}


#This is starting to overlap with the BaseImporter/Importer a little
#but new style Peak import is not supported by Importer just yet.

#todo expose generate_slices_from_names args
# change this to sub to use cache_slices, slices, slice_cache 
# there is some support for this in the BaseImporter,
# So we should probably move this code somewhere useable by all

#Fix this here for now, then move everything to the Helper, so we can access it
#here and from the Base/Importers

#get_Slice does not handle seem to filter skip_slices/slices either

sub slice_objects {
  my $self          = shift;
  my $slice_objects = shift;
  
  if($slice_objects){
    assert_refs($slice_objects, 'Bio::EnsEMBL::Slice', 'slice_objects');   
  }
  elsif(! defined $self->param_silent('slice_cache') ){
    $slice_objects = generate_slices_from_names
                      ($self->out_db->dnadb->get_SliceAdaptor, 
                       $self->slices, $self->skip_slices, 'toplevel', 
                       0, 0, $self->assembly); 
    #0, 0 are non_ref and inc_dups flags
    
    #Check we have some here
  }
   
  if($slice_objects){
    
    if($self->param_silent('slice_cache')){
      warn "Over-writing existing slice_objects";
      $self->param_silent('slice_cache', {});  
    }
    
    my %slices;
    
    map {$slices{$_->seq_region_name} = $_} @$slice_objects;

    $self->param('slice_cache', \%slices); 
  } 
   
  return [ values %{$self->param('slice_cache')} ]; 
}


#to do validate assembly, and that start = 1 if we get a slice name?

sub get_Slice {
  my ($self, $seq_region) = @_;

  #in case we have passed a slice name
  if($seq_region =~ /^\S*:\S*:(\S+):\S+:\S+:\S/) { $seq_region = $1; }
  #should we validate assembly here?
  
  #In case UCSC input is used... carefull names may not match with ensembl db!
  $seq_region =~ s/^chr//i;   #case insensitive     

  #We have seen a slice, but have not restricted slices so this
  #must be a slice we can't handle
  my $slice = undef;
  
  $self->slice_objects if ! defined $self->param_silent('slice_cache');
  my $slice_cache = $self->param_required('slice_cache');
  
  if( (! exists $slice_cache->{$seq_region}) &&
      (! ($self->slices || $self->skip_slices) ) ){
    
        
    #we're not actually testing slices and skip slices here
    #These shoudl ideally be cached, so we can access them quickly
    #or just ignore them?    
        
    throw("Unable to get Slice for:\tx".$seq_region.'x');      
  }
  else{
    $slice =  $slice_cache->{$seq_region}; 
  }
     
  return $slice;
}



#todo
#Add input here?
#There is a potential clash between DataSet and ResultSet support here

#This won't be recalled nicely
#and will be unclear in the caller whether all set methods will be injected
#although each analysis 'should' know which methods will be available
#due to the dataf flow input

sub fetch_Set_input{
  my ($self, $return_set_type) = @_;
  
  if(! defined $return_set_type){
    thow("Return set type argument is not defined, must be one of:\n\t".
      "DataSet FeatureSet ResultSet InputSet");
  }
  elsif($return_set_type !~ /^(Result|Feature|Data)Set$/o){
    thow("$return_set_type is not a valid return set type, must be one of:\n\t".
        "DataSet FeatureSet ResultSet InputSet");
  }  
  
  #why do we have a mistmach between the case of the set_type and the method name?
  #revert to title case i.e API standard as opposed to table name
  #as these will create methods!
  
  my $set_type       = $self->param_required('set_type');
  my $adaptor_method = 'get_'.$set_type.'Adaptor'; 
  my $dbid           = $self->param_required('dbID');
  my $db             = $self->param_required('out_db');
  my $set_name       = $self->param_required('set_name');
  #can't $db->can($adaptor_method) as this doesn't work with autoload
  
  
  $self->helper->debug(1, "Fetching $set_name $set_type with dbID $dbid");
  my $set            = $db->$adaptor_method->fetch_by_dbID($dbid);
  
  if(! defined $set){
    throw("Could not fetch $set_type with dbID $dbid ($set_name)"); 
  }
  elsif($set->name ne $set_name){
    throw("Fetch $set_type with dbID $dbid, expected $set_name but got ".$set->name);  
  }
  
  $self->set_param_method($set_type, $set);
  
  if($set_type eq 'DataSet'){
    my @rsets = @{$set->get_supporting_sets('result')};
    
    if(scalar @rsets != 1){
      
      if($return_set_type eq 'ResultSet'){
        throw("Expected 1 ResultSet, found:\n".$self->helper->dump(\@rsets));
      }
    }
    else{
      $self->helper->debug(2, "Setting result_set:\t".$rsets[0]->name);
      $self->set_param_method('ResultSet', $rsets[0]);
    }
  }
  elsif($set_type eq 'ResultSet'){
    my @dsets = @{$db->get_DataSetAdaptor->fetch_all_by_supporting_set($set)};
    
    if(scalar (@dsets) != 1){
      
      if($return_set_type eq 'DataSet'){
        throw("Failed to identify unique DataSet using ResultSet:\t".$set_name);
      } 
    }
    else{
      $self->set_param_method('DataSet', $dsets[0]);
      $self->helper->debug(2, "Setting data_set:\t".$dsets[0]->name);  
    }
  }
  else{
   throw($set_type.' set_type not supported. Must be DataSet, ResultSet or InputSet'); 
  }
  
  
  
  if($self->param_silent('DataSet')){
     my $fset = $self->DataSet->product_FeatureSet;
     
     if($fset){
        $self->set_param_method('FeatureSet', $fset); 
        $self->helper->debug(2, "Setting feature_set:\t".$fset->name);        
     }
  }
  
  if( ($return_set_type eq 'FeatureSet') &&
      (! $self->param_silent('FeatureSet')) ){
    throw("Failed to fetch a FeatureSet using $set_type:\t".$set_name);
  }
  
  #if we don't specify return_set_type
  #then an analysis may get an unexpected set returned if the data flow isn't correct
  #hence needs to be mandatory
  
  return $self->param_required($return_set_type);
}

=over

=item data_set

=item result_set

=item feature_set

=item slice_objects

Only if they have been specified or are required

=cut

#Do same for input_set_name(no, this will just be the set name) and result_set_name?



#todo move/remove the rest here

#These are specific to batch job
#Move to BaseBatch?
#Or move IdentifySetInput array methods to here a different base class (factory?)
#what else will be using feature_types etc. Need to see pipeline diagram!


#these are only needed when adding the initial input_ids
#analysis specific batch jobs can probably just use analysis


#There are basically two types of jobs
#One which uses arrays of meta data to identify sets

#Then serial or batch jobs which take one or more input_ids of dbID and name
#but these jobs will have access to the pipeline_wide_params from meta
#only redefine with generic names if the runnable can be generic also (e.g.IdentifySetInputs)


#These are always optional as they have defaults

#TO DO remove these as they will be replaced with get/set_param_method in Runnables
#where we can specify (ion context) whether they are required or not

#sub feature_set_analysis { return $_[0]->param_silent('feature_set_analysis'); }

#sub input_set_analysis   { return $_[0]->param_silent('input_set_analysis');   }

#This has been superceded by more specific alignment_analysis 
#we may have want to create different result_sets with different analyses
#this translation between params and result_sets shoudl be done by the runnables
#so we should really have generic set analyses here
#TO DO, rename feature_set_analysis peak_analysis!
#sub result_set_analysis  { return $_[0]->param_silent('result_set_analysis');  }


#Are these mandatory if accessed?
sub cell_type    { return $_[0]->param_silent('cell_type');} #optional!

sub feature_type { return $_[0]->param('feature_type'); }

sub experimental_group   { return $_[0]->param('experimental_group');           }






1;
