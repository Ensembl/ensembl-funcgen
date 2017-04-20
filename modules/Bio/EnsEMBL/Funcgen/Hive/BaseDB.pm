package Bio::EnsEMBL::Funcgen::Hive::BaseDB;

use warnings;
use strict;

use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor; 
use Bio::EnsEMBL::Funcgen::DBSQL::TrackingAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( generate_slices_from_names assert_refs );
use Bio::EnsEMBL::Utils::Exception         qw( throw warning );

use base ('Bio::EnsEMBL::Funcgen::Hive::Base');

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
  
#   $self->validate_dir_param
#    ('db_output_dir',
#     1,  #create/writeable? #should this be undef, and force the runnables to validate_dir_param output_dir?
#     $self->data_root_dir.'/output/'.$self->out_db->dbc->dbname); #default
  
  
  $self->process_params(['slices', 'skip_slices'], 1, 1);#optional/as array flags
 
  #validate/generate slices/skip_slices here
  if($self->slices || $self->skip_slices){
    $self->slice_objects; 
    #Validates custom slice & skip slices before we start fanning jobs
  }

  return;
}

sub slice_objects {
  my $self          = shift;
  my $slice_objects = shift;
  my $param_dups    = $self->param_silent('include_slice_duplicates');

  if($slice_objects){
    assert_refs($slice_objects, 'Bio::EnsEMBL::Slice', 'slice_objects');   
  }
  elsif(! defined $self->param_silent('slice_cache') ){
    $slice_objects = generate_slices_from_names
                      ($self->out_db->dnadb->get_SliceAdaptor, 
                       $self->slices, $self->skip_slices, 'toplevel', 
                       0, $param_dups, $self->assembly); 
    #0, 0 are non_ref and inc_dups flags
  }
   
  if($slice_objects){
    
    if($self->param_silent('slice_cache')){
      warn "Over-writing existing slice_objects";
      $self->param_silent('slice_cache',    []); 
      $self->param_silent('slice_registry', {});  
    }
    
    #Deref for safety, so we don't get any odd updating of this cache
    #via the passed arrayref.
    my $slices = [@$slice_objects];
    $self->param('slice_cache', $slices);
    #Point the registry at the cache to reduce footprint 
    my %slice_registry;
  
    foreach my $i (0..$#{$slices}){
      $slice_registry{$slices->[$i]->name} = $slices->[$i];

      my $sr_name = $slices->[$i]->seq_region_name;

      if(! exists $slice_registry{$sr_name}){
        $slice_registry{$sr_name} = $slices->[$i];
      }
      else{ #Create and array to handle PARs

        if(ref($slice_registry{$sr_name}) ne 'ARRAY'){
           $slice_registry{$sr_name} = [$slices->[$i]];
        }

        push @{$slice_registry{$sr_name}}, $slices->[$i];
      }
    }

    $self->param('slice_registry', \%slice_registry);
  } 
   
  return $self->param('slice_cache'); 
}




sub get_Slice {
  my ($self, $name, $lstart, $lend) = @_;

  $self->slice_objects if ! defined $self->param_silent('slice_registry');
  my $slice_cache = $self->param_required('slice_registry');

  #In case UCSC input is used... carefull names may not match with ensembl db!
  #This should not alter slice names
#   $name =~ s/^chr([^o])/$1/i;   

  if( (! exists $slice_cache->{$name}) &&
      (! ($self->slices || $self->skip_slices) ) ){    
    #We have seen a slice, but have not restricted slices so this
    #must be a slice we can't handle     
    throw("Unable to get Slice for:\t$name");      
  }

  my $slice = (exists $slice_cache->{$name}) ? 
   $slice_cache->{$name} : undef;

  if($slice && 
     (ref($slice) eq 'ARRAY')){
    #We have non-PAR regions cached as the slice cache was likely generated
    #without no_dups set.

    if(! ($lstart && $lend) && 
       ! ($self->param('include_slice_duplicates') )
      ){
      throw('The slice cache has been generated without duplications i.e. '.
        'it contains multiple Y non-PAR slices. This requires specifying a '.
        'loci start and end value to resolve the slice required.'.
        "\nAlternatively, generate the slice cache with the inc_dups argument set\n".
        "Note: Setting inc_dups will result in duplicate features fetched across the X-Y PAR regions.");
    }

    #Identify the appropriate non-PAR slice
    my $tmp_slice;

    foreach my $non_par(@{$slice}){

      if(( ($lstart >= $non_par->start) && ($lstart <= $non_par->end)) ||
         ( ($lend   >= $non_par->start) && ($lend   <= $non_par->end)) ){
        $tmp_slice = $non_par;
        last;
      }
    } 

    if(! defined $tmp_slice){
      throw("Failed to find a non-PAR slice for:\t$name $lstart - $lend\n".
        'This location must lie on a PAR region or outside the sequence');
      #Print non_PARs here?
    }

    $slice = $tmp_slice;
  }
     
  return $slice;
}

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
  
  # inject these also with set_param_method?  
  my $set_type       = $self->param_required('set_type'); 
  my $dbid           = $self->param_required('dbID');
  my $db             = $self->param_required('out_db');
  my $set_name       = $self->param_required('set_name');
  #can't $db->can($adaptor_method) as this doesn't work with autoload
  
  my $adaptor_method = 'get_'.$set_type.'Adaptor'; 
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

#Are these mandatory if accessed?
sub epigenome    { return $_[0]->param_silent('epigenome');} #optional!

sub feature_type { return $_[0]->param('feature_type'); }

sub experimental_group   { return $_[0]->param('experimental_group');           }






1;
