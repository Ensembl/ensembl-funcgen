
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

Bio::EnsEMBL::Funcgen::Hive::Base

=head1 DESCRIPTION

A base class for all other Funcgen hive modules. Provides general methods
for pipeline infrastructure, config processing and validation.

=cut

package Bio::EnsEMBL::Funcgen::Hive::Base;

use warnings;
use strict;

use Devel::Peek;

use Bio::EnsEMBL::Funcgen::Utils::Helper;
use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( scalars_to_objects 
                                               validate_path
                                               get_files_by_formats 
                                               path_to_namespace );
use Bio::EnsEMBL::Utils::Scalar            qw( assert_ref check_ref );  
use Scalar::Util                           qw( blessed );                                            
                                               

use base qw(Bio::EnsEMBL::Hive::Process);



# TODO
# Move a lot of this to BaseSequenceAnalysis
# what won't be using the tracking DB? i.e. Isn't everything a BaseDB?
#We can do peak calling and alignment outside of DB using -use_tracking_db 0

#This means we need to separate the requirement for tracking DB from BaseDB as make it available
#via Base too

# How can jobs which don't use any DB (tracking or output) get access to the necessary params?
#simply they can't? local file locations are stored in the tracking tables
#No we can allow manual definition of input_ids containing the file paths, hence bypassing IdentifySetInputs
#Can we even do this with the motif feature mapping too? Or does that depend on the DB? Could probably make it independant
#maybe by dumping data from the DB and using that as input. Hence users could simply recapitulate the dumps as a work around
#and define some input parameter to point to them
#This would change the structure of the pipeline somewhat, could probably be done with a different config (skipping IdentifySetInput and DefineSets)
#Need to be mindful when developing runnables
#basically all runnables DB interactions would have to be depdendant on use_tracking_db
#when we normally load into the DB, we would need to check for the db param or some -no_output_db flag
#this might be useful if we want to run an analysis but don't want it to write to the DB for some reason

# Go ahead with use_tracking_db = 0 implmentation, 
# and code runnables to support no_output_db implementation, i.e. take DB or non DB input params 
# but leave the non db config for now. Infact first analysis is StoreRollbackSets, so this would have to be skipped 


#global values for the Helper... maybe pass as parameters...
#$main::_debug_level = 0; Now set below
$main::_tee         = 0;
$main::_no_log      = 1;


#Used in set_param_arrays for scalars_to_objects
#Could our this to make it dynamically mutable, i.e. redefine it from a sub class
my %param_class_info = 
  (
    cell_types          => ['CellType',          'fetch_by_name'],
    feature_types       => ['FeatureType',       'fetch_by_name'],
    analyses            => ['Analysis',          'fetch_by_logic_name'],
    experimental_groups => ['ExperimentalGroup', 'fetch_by_name'],
    cell_type           => ['CellType',          'fetch_by_name'],
    feature_type        => ['FeatureType',       'fetch_by_name'],
    analysis            => ['Analysis',          'fetch_by_logic_name'],
    experimental_group  => ['ExperimentalGroup', 'fetch_by_name'],
  );


my %object_dataflow_methods = 
  (
   'Bio::EnsEMBL::Analysis' => 'logic_name', 
  );



my %valid_file_formats = 
 (
  bed  => 'Bed',
  sam  => 'SAM',
  bam  => 'BAM',
 );

#Advantage of having separate fetch_input, run and write methods is to allow
#calling of super methods at appropriate point.
#Theoretically, this could all be done in the run method with calls to the super methods at the 
#appropriate points, but this would break the object model and make it unihneritable.

#This defines a set of parameters based on given parameters from the pipeline:
sub fetch_input {   # nothing to fetch... just the DB parameters...
  my $self = $_[0];
  $self->SUPER::fetch_input;
  
  $main::_debug_level = $self->debug; #Pass this onto the Helper and any other modules in the process

  #Catch and validate generic params here
  #but manage mandatory aspect in more specific runnables
  

  
  #validate_dir_params uses process_params which uses get_param_and_method
  #Hence all these dir with now have a method available
  
  #root_input_dir is mandatory and workdir and hive_output_dir have defaults based on this
  #Mandatory
  $self->validate_dir_param('data_root_dir'); 
 
  $self->validate_dir_param('bin_dir',    undef, 1); #optional

  #Optional with defaults
  $self->validate_dir_param
   ('hive_output_dir',    
    1, #create/writeable
    $self->data_root_dir.'/output/'.$self->param_required('pipeline_name').'/hive_debug'
   ); 
    
  $self->validate_dir_param
   ('work_root_dir', 
    undef, #Let's not create this until we need it
    $self->data_root_dir.'/output/'.$self->param_required('pipeline_name')
   );
  
  #output and work dir need setting based on dbname
  if($self->get_param_method('use_tracking_db')){
    $self->_set_out_db;      
  }

  return;
}



#These are already defined in Hive::Process
#sub write_output { my $self = shift @_;  return 1; }
#sub run          { my $self = shift @_;  return 1; }

# Validate mandatory non DB params separately 
# as they are not mandatory for subclass BaseDB
# This should be called for any cubclass which is not BaseDB

sub validate_non_DB_inputs{
 my $self = $_[0];
 
 $self->get_param_method('species', 'required');
 $self->get_param_method('assembly', 'required'); 
 return;
}


#Getter methods to catch potential typos at compile time rather than runtime or not at all!
#Use param directly as the previous set method was not catching undef (as param only warns with undef)


sub cell_type { return $_->param_required('cell_type'); }


#add other generic params here
#cell_type feature_type etc


#This does allow config over-ride but would need adding to the pipeline_wide_params

sub alignment_root_dir {
 my $self = $_[0];
 
 if(! $self->param_silent('alignment_root_dir')){
   $self->set_dir_param_method('alignment_root_dir', [$self->data_root_dir,
                                               'alignments',
                                               lc($self->param('species')),
                                               $self->param_required('assembly')]);
                                               
   #complete path will include study/experiment name and input_set logic_name
   #which will be the logic name of the alignment
 } 
  
 return $self->param('alignment_root_dir'); 
}

#todo support >1 root_output dir?
#might want to spread data across two scratch areas?
#move alignments to output dir?
#change to just root_data_dir
#no need for root_output_dir really as this only supports output_dir (e.g. output/dbname)

#This will be required by the alignments and the collections
#Currently hardcodes the alignment dir based on the alignment_root_dir
#Takes input_set as this is the return type of ResultSet::get_support


#create no longer required in here as we normally 
#create this using get_output_work_dir which does the validation

#TODO add support for control alignment dir!

sub alignment_dir {
  my ($self, $rset, $create, $control) = @_;
   
  if(defined $rset){ 
    #my $exp_id = $rset->get_Experiment->dbID;
    
    #todo handle input_set analyses here
    #as the logic_name will eventually be added to the path
    #hence we can't have mixed analyses
    #for my $i(1..$#$input_sets){
    #  if($exp_id != $input_sets->[$i]->get_Experiment->dbID ){
    #    throw('Cannot define a single alignment_dir for InputSets which relate '.
    #          "to different Experiments:\n".$self->helper->dump($input_sets)); 
    #  }
    #}
    
    
    #This is not quite the same format as the old alignment_dir
    #we are putting all InputSets in the same location
    #which is fine
    
    $self->set_dir_param_method('alignment_dir', 
                                [$self->alignment_root_dir, 
                                 $self->get_study_name_from_Set($rset, $control)],
                                $create);    
  }
  return $self->param('alignment_dir');
}



sub get_study_name_from_Set {
  my ($self, $set, $control) = @_;
  
  if(! (ref($set) && 
        ($set->isa('Bio::EnsEMBL::Funcgen::InputSubset') ||
         $set->isa('Bio::EnsEMBL::Funcgen::ResultSet') ))){
    throw('Must pass a valid InputSet or ResultSet');         
  }
  
  my ($exp_name, $ftype);
  my $ctype = $set->cell_type->name;
  
  if($control){
    $control = _get_control_InputSubset($set);
    
    #This is based on the assumption that all non-ctrl subsets are associated
    #with one ResultSet/Experiment.
    
    my $exp   = $control->get_Experiment;
    $exp_name = $exp->name;
 
    foreach my $isset(@{$exp->get_InputSubsets}){
      
      if(! $isset->is_control){
        $ftype = $isset->feature_type->name;  
        last;
      }  
    }      
  }
  else{  
    $exp_name = $set->get_Experiment->name ||
      throw("Cannot find unique experiment name for ResultSet:\t".$set->name);
    $ftype = $set->feature_type->name;
  }
  
  (my $study_name = $exp_name) =~ s/${ctype}_${ftype}_(.*)/$1/;
  return $study_name;
}


#This seems odd as the naming convention is
#$experiment->name
#However, we need to hande controls here which will not match the experiment name
#and experiment name is also unsafe as experiment is really a study, and 
#need not contain the cell_type and feature_type
#Given that we are already making many assumptions based on this we can do the same here?
#although it would be good to have a set of methods to handle this
#such that we don't pepper the API with assumptions i.e. if we want to change it
#it only changes in one location.

#So we need a method which will generate the standard set prefix and the control set prefix
#based on an existing Set, or an Experiment. 
#Can't d this for Experiment with an additional cell_type/feature_type being passed
#which means that the implementations would still have to change should we update the 
#nomenclature


#Should probably pass prefix name generation off to another method 
#taking all the string elements. Then we would have one method which 
#defines the nomenclature, and would not be dependant on API object
#Would still require rejigging of implementation should the nomenclature
#change, and API object are advantageous here for validation
  

#Don't bother validating that all the controls have matching ctypes 
#and ftyps as this should be done during import
#Just assume we can pick 1 which is representative


#Will fail if we allow inter group controls?
#No, we should never have mixed controls
#Check this is the case on import
    
sub _get_control_InputSubset{
  my $set = shift; 
  my @is_sets;
  
  if( $set->isa('Bio::EnsEMBL::Funcgen::ResultSet') ){
    @is_sets = @{$set->get_support('input_subset')};
  
    if(! @is_sets){
      throw("Failed to identify control InputSubset support for ResultSet:\t".$set->name);  
    }
  }
  else{
    @is_sets = ($set);
  }
  
  my @ctrls = map { $_ if $_->is_control } @is_sets;

  if(! @ctrls){
    throw('Could not identify a control InputSubset from '.ref($set).":\t".$set->name);  
  }

  return $ctrls[0];
}
  

sub get_set_prefix_from_Set{
  my ($self, $set, $control) = @_;
  
  my $study_name = $self->get_study_name_from_Set($set, $control);
 
  if(! (ref($set) && 
        ($set->isa('Bio::EnsEMBL::Funcgen::InputSubset') ||
         $set->isa('Bio::EnsEMBL::Funcgen::ResultSet') ))){
    throw('Must pass a valid InputSet or ResultSet');         
  }
  
  my $ftype;
  
  if($control){
    $ftype = _get_control_InputSubset($set)->feature_type->name;
  }
  else{
     $ftype = $set->feature_type->name;
  }
 
  return $set->cell_type->name.'_'.$ftype.'_'.$study_name; 
}


sub get_output_work_dir_methods{
  my ($self, $default_odir, $no_work_dir) = @_;
   
  my $out_dir = $self->validate_dir_param('output_dir', 1, $default_odir); 
  my $work_dir;    
             
  if(! $no_work_dir){
    my $dr_dir             = $self->data_root_dir;
    
    if($default_odir !~ /$dr_dir/){
      throw('Cannot set a work_dir from an output dir('.$out_dir.
        ") which is not in the data_root_dir:\n\t$dr_dir");   
    }    
               
    my $wr_dir             = $self->work_root_dir;
    ($work_dir = $out_dir) =~ s/$dr_dir/$wr_dir/;
  
    if($work_dir eq $out_dir){
      $work_dir = $wr_dir;
      warn("Failed to create appropriate work_dir subdirectory, defaulting to work_root_dir:\n\t".
        $work_dir);
    }
    
    $self->validate_dir_param('work_dir', 1, $work_dir);
  }
  
  $self->helper->debug(1, "output_dir $out_dir\nwork_dir $work_dir\ndefault $default_odir");

  return ($out_dir, $work_dir);
}



#This is for useful if we want to output logs in run
#instead of write to DB in write
sub no_write     { return $_[0]->param_silent('no_write');  } 



#Although this is not BaseDB we can still use a tracking DB
#for analysis which don't have an explicit (non-tracking) output to the DB e.g. alignments
#was originally called db, but looks like this redefined a Process method
#although there was no error output for this


#Todo move this to the importer also?
#This would mean all DB dependant runnables would need an importer
#but the importer currently requires a parser/vendor
#which is not always appropriate


sub _set_out_db {
  my $self = $_[0];
  
  my $db = $self->get_param_method('out_db', 'required');
  #This will either be the params hash or a DBADaptor
 
  
  #This isn't catching out_db being ill-defined in the config. but this would 
  #require extra control flow for every
 
  if(! ref($db)){ #Has been set to scalar in config?
     throw("out_db config is not a Hashref of DBAdaptor::new parameters:\t$db");
  }  
  else{
    
    if(ref($db) eq 'HASH'){
      my $dnadb_params = $self->param_required('dnadb');  
      #Always required as we don't want to default to ensembldb for a pipeline!
          
      #Create TrackingAdaptor here, as we can't get_TrackingAdaptor later
      my $adaptor_class = ($self->use_tracking_db) ? 
        'Bio::EnsEMBL::Funcgen::DBSQL::TrackingAdaptor' :
        'Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor';
       
      eval{	$db = $adaptor_class->new(%{ $db }, %{ $dnadb_params }) };	  
      if($@) { throw "Error creating the Funcgen DBAdaptor and/or dna DBAdaptor $@";  }  
      
      #Reset the $db to an actual DBAdaptor (rather than a Tracking/BaseAdaptor)
      if ($self->use_tracking_db){
        $self->set_param_method('tracking_adaptor', $db);
        $db = $db->db; 
      }
        
  	
      #Actually test connections
      $db->dbc->db_handle;
      $db->dnadb->dbc->db_handle;
  	
      #To avoid farm issues...
      if($self->param('disconnect_when_inactive')){
        $db->dbc->disconnect_when_inactive(1);
        $db->dnadb->dbc->disconnect_when_inactive(1);
      }
      
          
      # VALIDATE/SET Assembly
      #This may clash with default assembly from DB if we are loading onto an old assembly
      #normally we would project this across during import, but it might be valid to 
      #actually load on an old assembly for comparison
      #This can't be overloaded for assembly and new assembly as this will clash in the original input ids  
      #Parameter should be -projection_assembly or -target_assembly?
      my $assembly   = $self->get_param_method('assembly', 'silent');
      my $cs_adaptor = $db->dnadb->get_CoordSystemAdaptor;
    
      if(! defined $assembly){
        $assembly = $cs_adaptor->fetch_by_rank(1)->version; 
        
        if(! $assembly){
          throw("Failed to identify default assembly from the DB");
        }
        
        $self->assembly($assembly);
      }
      else{ #validate it exists in the dnadb
         
        my $cs = $cs_adaptor->fetch_all_by_version($assembly)->[0];
      
        if(! defined $cs){
          throw("The -assembly $assembly does not exist in the database:\t".
            $db->dnadb->dbc->dbname); 
        }
        #todo make sure all Slice/CoordSystem calls use this assembly explicitly.
      }
     
      # VALIDATE/SET Species
      my $param_species = $self->get_param_method('species', 'silent');
      my $db_species = $db->species;
      
      if( ($param_species && $db_species) &&
        (lc($param_species) ne lc($db_species)) ){
          throw('Mismatch between the DB species (meta species.production_name) and'.
          " -species parameter:\t".$db_species."\t".$param_species);      
      }
      elsif(! defined $param_species){
      
        if(! defined $db_species){
          throw('Must either set -species or use a DB where the correct species '.
            'name is set in the meta table as species.production name');  
        }
      
        # Have to be mindful that we can't access the adaptors from the 
        #Registry using the -species parameter. Hence it is unsafe to use the Registry!
        #Always get the adaptors directly from the DB.
        $self->species($db_species);
      } 
      #else we don't care so much about it not being set as db species/species.production name?
      #if it isn't then registry based access won't work, as we know.      
    }
    elsif(! check_ref($db, 'Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor')){
      throw("The out_db param is set to an unexpected reference:\t$db\n".
        "Please check the out_db config, which should be a Hashref of DBAdaptor parameters");
    } 
  }
  
  return $self->out_db($db);
}



#Expose param_silent for handling optional params

sub param_silent {
  my $self = shift;
  return $self->input_job->_param_silent(@_); 
}


#Move this to Base? Or just merge this whole module with base.
#change this to set_param
#this will handle scalars and arraysref to objects or arrays/arrays of objects
#This could also handle hashref which are a constraints_hash to pass the the relevant
#fetch_all method. This would require the requisite support from fetch_all
#so we have to catch that here, as there is no way of detecting that.

#todo add skip fetch and alt class arrays ref args?


sub process_params {
  my ($self, $param_names, $optional, $as_array) = @_;
  
  #handle single scalar param name
  my $ref_type = ref($param_names);
  if($ref_type){ #can be ''
    
    if($ref_type ne 'ARRAY'){
      throw("Param names argument cannot be a $ref_type.\n".
        'Only a single scalar or an Arrayref of param names is permitted.');  
    }
    
    if(! scalar @$param_names){
      throw("Cannot pass an empty Arrayref for param names argument\n".
        'Only a single scalar or an Arrayref of param names is permitted.'); 
    }
    
  }#Don't need to catch undef or null string as this will be caught by param method
  else{ #Must be single true scalar param name
    $param_names = [$param_names] 
  }
   
 
  my $param_method = ($optional) ? 'silent' : 'required';
  my %all_params;
  
  foreach my $param_name(@$param_names) {
    my $param_val = $self->get_param_method($param_name, $param_method);
   
    #defined here will not catch empty hashes
    #arrays are not permitted here
    
    
    if(defined $param_val){
      my $ref_type = ref($param_val);
      my $values = $param_val;
      
      if(! $ref_type){
        if($as_array){
          $values = [split(/,/, $values)];
        }
      }
      elsif( ($ref_type ne 'HASH') &&
              ($ref_type ne 'ARRAY') ){
        throw('process_params can only handle scalar (comma separated), '.
          'Arrayref of scalars or Hashref values');  
      } 
      
      #warn "before param_class_info test with $param_name";
             
      if(exists $param_class_info{$param_name}){
        warn "yay $param_name exists in param_class_info";
 
        if($ref_type eq 'HASH'){#can only be hash
          throw("Need to implement $param_name fetch_all for process_params");
        }
        else{ #We have an Arrayref of scalars
          $values = scalars_to_objects($self->out_db, @{$param_class_info{$param_name}},
                                       $values);
        }
        
        $self->helper->debug(1, "$param_name is now ", $values);
      }
      
      
      #warn "Setting $param_name to $values";
      
      $all_params{$param_name} = $values;
      $self->$param_name($values); 
    }
    else{ #can only be optional $optional
      #This may already be set to undef, but just in case set here
      #so we don't keep triggering the warnings!
      #$self->param($param_name, undef); #This will not work!
      $self->input_job->{'_param_hash'}{$param_name} = undef;
      
      #This won't cause issues if we later want to param_required as it tests if defined
      #Do we want to set this in all_params?
    }
  }
  
  return \%all_params; 
}



sub validate_dir_param{
  my($self, $dir_name, $create, $optional_or_default) = @_;
 
  my $req_or_silent = 'required';
  my $default;
 
  if($optional_or_default){
    $req_or_silent = 'silent';   
  
    if($optional_or_default ne 1){
      $default = $optional_or_default;  
    }  
  }
 
  
  my $path = $self->get_param_method($dir_name, $req_or_silent, $default); 
   
  #This will have throw if is it required and not defined or has default
 
  validate_path($path, $create, 1, $dir_name);
 
  #my %dirs = %{$self->process_params($dir_param_names, $optional)};
  #foreach my $dir_name(keys %dirs){ 
  # my $path = $dirs{$dir_name};

  #  if(defined $path){
  #    validate_path($path, $create, 1, $dir_name);  #1 is dir flag
  #  }
  #}
   
  return $path;  
}

#$path can be array ref to pass to FileSpec or just a string

#rename this to set_dir_param_method

#separate this?
#currently calling this from pre-exiting dir methods to set the default dir
#see alignment_dir

sub set_dir_param_method{
  my ($self, $dir_param_name, $path, $create) = @_;      
  
  if(! (defined $dir_param_name &&
        defined $path) ){
    throw('You must provide directory parameter name and path arguments to validate');
  }
  
  
  #$self->helper->debug(1, "set_dir_param $path, $create, 1, $dir_param_name", $path);
  
  $path = validate_path($path, $create, 1, $dir_param_name);#1 is dir flag
  
  $self->helper->debug(1, "Setting $dir_param_name:\t$path");
  
  return $self->set_param_method($dir_param_name, $path);
}


sub helper {
  my $self = $_[0];
  #allow overwriting of $main::_tee/no_log/_debug_level?
 
  my $helper = $self->param_silent('helper');

  if(! defined $helper){
    $self->param('helper', Bio::EnsEMBL::Funcgen::Utils::Helper->new());
  }
  
  #todo CR enable help param defs e.g. debug level etc
  
  return $self->param('helper');
}




=head2 get_param_method

  Arg [1]    : String - param name to generate method for.
  Arg [2]    : String - Optional param method type: 'silent' or 'required'
               Default will just use 'param'.
  Example    : $self->get_param_method('some_parameter');
               my $param = $self->some_parameter; 
  Description: Defines a wrapper method to 'param' for the given parameter name, 
               and returns the parameter via the defined param method type.
               This prevents usage of strings when using the param method
               and the danger of uncaught typos. It also protects against
               redefining pre-existing methods by checking if they already 
               exist. This is very useful as debugging can be painful when you
               have trashed a hive API method unknowingly
               WARNING: This does obfuscate the codebase as dedicated methods 
               are removed. This can be ameliorated by the addition of POD for 
               those methods which are injected by get_param_method. <-- DO THIS!!!
  Returntype : Dependant on the config/input_id.
  Exceptions : Throws if method is already defined
               Throws if param method type is not valid
  Caller     : General
  Status     : At risk

=cut


#Issues
#It is still fairly hard to know when a method may have been injected
#this maybe conditional on subclass calling a Base method etc.

#todo utilise this in process_param methods?

#could even put in a sneaky grep to ensure that the pod has been added?

#This doesn't account for a param being required after it is first set

#todo use AUTOLOAD to catch methods which may not have been injected
#and give a useful warning?

#todo add more wrappers for silent and required, such that we don't need to 
#pass string args, which will catch errors at compile rather than runtime
#args are fully validated so this isn't unsafe

#can't have default if required???
#would throw in param_required if not defined 

sub get_param_method {
  my ($self, $param_name, $req_or_silent, $default) = @_;
  my $value = $self->_param_and_method($param_name, undef, $req_or_silent);
  
  if((! defined $value) && 
     defined $default ){
    #warn "There is not $param_name defined in the config. Defaulting to:\t$default";   
    $self->$param_name($default);
    $value = $default;     
  }
   
  return $value;  
}

sub set_param_method {
  my ($self, $param_name, $param_value, $req) = @_;

  if(defined $req){
    
    if($req ne 'required'){
      #do this instead of a boolean, to avoid someone specifying silent and getting unexpected
      #behaviour, also mirrors get_param_and_method interface
      throw("$req is not a valid method type for set_param_and method, can only be 'required'");  
    }
    elsif(! defined $param_value){
      throw($param_name.' value is required but not defined');  
    }
  }
  
  return $_[0]->_param_and_method($param_name, $param_value, $req);    
}


#This is basically Class::Accessor but tied into the hive system
#and validates we are not over-writing some other method
#todo improve pacakge_method test? Should we ever need to over-write and
#existing method?

sub _param_and_method {
  my ($self, $param_name, $param_value, $req_or_silent) = @_;
  my $cref;
  
  #First check the method does not already exist
  if($cref = $self->can($param_name) ){  
    
    my $package_method = &CvGV_name_or_bust($cref);
    
    #warn "CvGV name $package_method";

    if(defined $package_method){
      #Should not redefine non-Funcgen method
  
      if($package_method !~ /^Bio::EnsEMBL::Funcgen::/){
        throw("Caught attempt to redefine non-Funcgen method:\t$package_method\n".
          "Please use an alternative method name, or use the above method directly if appropriate");
      }
    }
    #else it does exist and we are likely simply calling this
    #method in a batch job, where the global scope has persisted from
    #the first job which initialised the method
    #add debug statement here?
    
    #We are still injecting/over-writing the coderef here
    #can we simply return is the package_method matches what we expect?
    #and throw if is doesn't
    
  }
  #else it does not exist, so safe to inject
    
   
  
  my $param_method = 'param';
  
  
  if(! defined $param_value){
  
    if(defined $req_or_silent){
    
      if( ($req_or_silent ne 'silent' ) &&
          ($req_or_silent ne 'required') ){
        throw("Param method type can only be 'silent' or 'required', not:\t$req_or_silent");    
      }
    
      $param_method .= '_'.$req_or_silent;
    }
  }
  #else we have already dealt with 'required' in set_param_and_method
  #'silent' is NA
  
  
  #This needs to use param_method
  #other wise we will get param warnings
  #if use silent and the param is not defined
  #and then call the injected method
  #or should we just use param_silent
  #and rely on the $pram_method call below
  
  
  if(! defined $cref){
    #Inject method, hopefully not Pragha Khan style!
    no strict 'refs';
  
    *{ref($self)."::${param_name}"} = sub { 
      my ($self, $param) = @_;
    
      if(defined $param){
        $self->param($param_name, $param);
      }
    
      return $self->param_silent($param_name);
    }; 
   
    use strict;
    #avoids being able to work with symbolic reference
    #i.e. Your::Package::Name->$param_method
    #rather than $obj->$param_method
  }
  
  #Use $param_method here rather than the $param_name method
  #we have just defined, so we can validate just once here
  #rather than every time the new method is called
  #This would cause problems if we reuse the method in another instance
  #within the same global scope(i.e. next batch job), 
  #as $req_or_silent may have changed
  
  
  #Can't pass undef $param value here as it will reset a defined value to undef!
  my @param_args = ($param_name);
  push @param_args, $param_value if defined $param_value;
  #to fx this properly and allow params to be set as undef, we would need to do the same trick
  #as Params::_param_silent, and use the size of @_ to see whether we truly had an undef
  #value passed or not
  #for now, will need to set undef directly using param
  
  
  return $self->$param_method(@param_args);
} 
  

  
  

# Here follows a list of param methods injected via set/get_param_and_method
# These are simple wrappers to the 'param' method 


=head2 get/set_(dir)_param__method defined 'param' accessor/setter methods

=over

=item data_root_dir

=item hive_output_dir

=item bin_dir

=item work_dir

=item use_tracking_db

=item out_db

Separated from _set_out_db so we don't have to do all the ref checking every time we access the db
Was just a getter, now is also a setter due to generic method implementation of get_param_and_method

=item species   

Defined by _set_out_db (if using a tracking DB or this is also a BaseDB) or validate_non_DB_inputs if not.

=item assembly

Defined by _set_out_db (if using a tracking DB or this is also a BaseDB) or validate_non_DB_inputs if not.

=back

=cut

#we definitely need branch_key_param|method
#for a truly generic implementation, i.e. we want to reuse module for a different analysis
#and need to branch based on a different method/param
#can we validate that here?
#or just do this in the caller?
#most params will have a dedicated method
#can do this whole validate by get_param_method in the caller.
#this will not validate method exists.
#we need to validate it only if branch_config has been defined

#branch_key_method still won't be enough if the return values should change as
#we iterate over values to data flow.
#will provide wrapper for nested method calls (as opposed to direct access)
#i.e. $self->$branch_key_method wraps $self->feature_set->analysis_logic_name
#This will prevent hardcoding of the branch_key, such that it can be redefined in the config.
#but will also require another wrapper if we don't have an existing direct access method

#sub init_branch_config {

#We should probably separate out the init_branching_by_method


#have branch_config wrapper method, which optionally takes ('analysis') or ('method', $method_name)
#as args, this then optionally calls the relevant init method?
#This is only useful if we want to get the config, when we are initing
#which we never really want?
#Too abstract, and actually obfucates initialisation


sub init_branching_by_analysis{  
  my $self = shift;
   
  #my $branch_config = $self->get_param_method('branch_config', 'silent');  
  #Not a passed param anymore, as we get it from the dataflow rules
 
 
  if(! defined $self->{branch_config}){
    
    my $dfr_adaptor = $self->db->get_DataflowRuleAdaptor;  
    
    if(! $dfr_adaptor->can('fetch_all_by_analysis_id')){
      #inject method here   
      no strict 'refs';
  
      *{ref($dfr_adaptor)."::fetch_all_by_analysis_id}"} = sub { 
        my $self    = shift;
        my $anal_id = shift;
        throw('Must provide and analysis id argument') if ! defined $anal_id;
        return $self->fetch_all("from_analysis_id=${anal_id}");
      }; 
   
      use strict;
    }
    
    if(! $dfr_adaptor->can('get_dataflow_config_by_analysis_id')){
      #inject method here using fetch_all_by_analysis_id 
      no strict 'refs';
  
      *{ref($dfr_adaptor)."::get_dataflow_config_by_analysis_id}"} = sub { 
        my $self    = shift;
        my $anal_id = shift;
        throw('Must provide and analysis id argument') if ! defined $anal_id;
        
        my %df_config;
        
        
        foreach my $dfr(@{$self->fetch_all_by_analysis_id($anal_id)}){
          #$anal_id here always represents the from analysis
          my $to_analysis = $dfr->to_analysis->logic_name;
          
          
          #Is it valid to wire to the same analysis using two different branches
          #We need to catch this and throw
          
          if(exists $df_config{$to_analysis}){
            throw('It appears that the pipeline configuration for '.$to_analysis.
              " has been wired via two separate branches:\t".$dfr->branch_code.
              ' & '.$df_config{$to_analysis}{branch}.
              "\nDynamic branch dataflow currently only supports 1 assoicated branch");  
          }
          
          $df_config{$to_analysis} ||= {branch => $dfr->branch_code, funnel=>undef };
          
          if($dfr->funnel_dataflow_rule_id){
            $df_config{$to_analysis}{funnel} = 
              $self->fetch_by_dbID($dfr->funnel_dataflow_rule_id)->to_analysis->logic_name;
          }
          
          #This will result in redundant branche entries across the to_analysis values
          #which is fine, we just need to handle this in the caller
          
          #How are we going to handle the potential of dataflowing down the same branch twice, due 
          #to the redundancy of the branches wrt to_analyis keys
          #this will warn, but not fail, and dataflow will not happen
          #so is safe-ish
          
          #This should be skipped over somehow
          #we could do this be caching refs to identify unique outputids
          #although there is nothing to stop the same output id being passed in a different array/hash
          #and hence a different ref
          #The way the runnable would and should be naturally written will prevent this
          #It will just be that a single analysis name will be used to flow to all?
          #we could add support branch_job_group by taking an array of analysis names?
          #so the args would be
          #$self->branch_job_group([[anal1, ...], [jobid1, ...], {funnel_analysis_name=>[funnel_job1, ...]}]);
          
          
          #branch_job_group will then check that all the analyses are on the same branch
          #and that they constitute all of the analyses i.e.
          #you don't add config without adding support and vice versa
          #this does not prevent the dynamic branch dataflow, just what is 
          #dataflown
          
          #Having the runnable directly linked to the logic names like this is 
          #normally a bad idea, but this is the sacriface made to enable
          #dynamic branch dataflow, and is better than having to proliferate
          #hardcoded dataflow code
          
          
          #dataflow_branched_job_groups b
        
        }
        
        return \%df_config;
      }; 
   
      use strict;
    }
    
    $self->{branch_config} = $dfr_adaptor->get_dataflow_config_by_analysis_id($self->analysis->dbID);
    
  }

  #if($validate_bkey_method){
  #  my $bkey_method = $self->get_param_method('branch_key_method', 'required');
    
      
    #warn "branch_key_method is $bkey_method";  
      
    #This will use an existing method or create one from en exiting parameter
    
  #  if(! defined $self->get_param_method($bkey_method, 'required')){
  #    throw("Required branch_key_method $bkey_method is not defined defined");
  #  }
    
    #Now reset branch_key_method method to actual branch_key_method?
    #or just go down the old 'method as a variable' root  
  #}
  
  #$branch_config ||= {};
  
  #100 is reserved branch for custom analyses?
  #Bin this af, as we may want to use 100 for something else
  #force config specification
  #my %branch_data_flow_ids = ();
  #don't set 1 and 100 now so we don't have
  #keys if we really don't flow anything
  #do we even need this now, as we have job groups
  #map { $branch_data_flow_ids{$_} = [] } values %$branch_config;
  
  #$self->helper->debug(1, "Initialised branch_data_flow_ids:", \%branch_data_flow_ids);
  
  #$self->set_param_method('branch_output_ids', \%branch_data_flow_ids);
    
  return $self->{branch_config};
}


#branch key is now the analysis name, and is part of the job group
#This does lose some of the generic nature of the branch config
#Will this work wrt IdentifySetInputs
#which is currently the only analysis whihc uses this?

#If branch is not defined, then job group contain branch_keys as analysis names
#If branch is defined, then job_group is simply and array of job_ids, with no branh keys
#

#$self->branch_job_group([{'Preprocess_BWA_replicate'=> [$job_id1, ...]}, {'RunIDR' => [$job_id2, ...]} ]);

#This method should validate that each analysis is either of the sme branch, or is a funnel for that analysis
#funnels have to come last? Or can we handle that?
#Should we convert analysis names to branch numbers here?


#Should we allow branching based on branch number only, when using branch config?
#We shouldn't have to conditionally call the branch method
#based on the analysis we are running, i.e. we should always dataflow
#The wiring or lack there of, should take care of this

#If we don't allow branching by number, then will we force use of analysis names
#which may appropriate for a given analysis (and they will not be wired)
#This will also create failures, if we are validating the analysis names
#

#We should still allow for branch_conf as this enable truly generic dataflowing
#without hardcoded analysis names


#Do not allow mixed branch number/analysis name job_grouping!
#This means we should never init_config, if we just going to branch by number
#so we should rename init_branch_config to init_branching_by_analysis!
#



#Delaying dataflow may be advantageous as it will prevent job flow until all 
#processing is done. In the case of creating a set for futher downstream analysis,
#the successful jobs may start flowing down stream. This then make rerunning the initial job
#problematic, as some of the inputs will have flown and some won't. Making setting of a
#rollback or recover mode problematic, as it will start rolling back a job which is has been successfully 
#flown. The only option here is to then figure out which jobs succeeded in flowing, and removing them 
#from the input, before resubmitting after correcting what caused the initial failure.
#If we delay the dataflow, then we can simply fix the problem, and retry the original job
#This is at the cost of delaying flow of jobs which don't have issues.
#The later is easier here.

#so we have two possible types of jobs group here
#one where the keys are analysis names
#and one where the keys are just branch numbers

#what about branch_key method jobs?
#the branch key method should be called in the runnable
#to pass the branck key here
#the init_branching_by_config method will have already set up
#the config. Does this allow funneling?
#Remember we aren't going to allow mixed branching
#or should we, so long as the key is not already defined in the config?

#let's deal with brancing by analysis/config first, which is effectively the same
#apart from the branching by config will not have any funnels.

#$self->branch_job_group([{'Preprocess_BWA_replicate'=> [$job_id1, ...]}, {'RunIDR' => [$job_id2, ...]} ]);
#This is currently not supporting multiple analyses for the same branch would have to change job_group format to following:
#$self->branch_job_group(['Branch2_analysis1', 'Branch2_analysis2'] [$branch2_job_id1, ...], {'RunIDR' => [$job_id2, ...]} ]);
#will only ever have 1 semaphore
#However, the reason why we are doing this dynamic braching (analysis specific resource managment) means we are very unlikely to 
#do this rather than just configuring another branch, unless they both feed into the same semphored job!
#no no no, this will not break the resource management as this is based on the target analysis, not the branch
#This will be to support two differeing analyses which require the same inputs
#we don't yet have an example of this, but it is valid


sub _get_branch_number{
  my $branch_codes  = shift;
  my $branch_config = shift;  
  assert_ref($branch_codes, 'ARRAY', 'Branch code array');
  my $branch;
 
  if(defined $branch_config){
    
    foreach my $bcode(@$branch_codes){
      
      if(! exists $branch_config->{$bcode}){
        throw("Could not find branch config for analysis or branch key:\t".$bcode.
          "\nDefault to a 'custom' branch here");
        #We would need and custom_analysis_name method
        #as we can't just use 'custom' as this may clash
        #and we don't know the analysis name template here
        #Could only support 1 custom analysis per runnable
        #unless we had a get_custom_analysis_name method
        #which would parse the unknown analysis name appropriately
        #overkill. stop.
      }
      else{
        $branch ||= $branch_config->{$bcode}{branch};
        
        if($branch_config->{$bcode}{branch} != $branch){
          throw($bcode.' fan analysis/branch key has branch '.$branch_config->{$bcode}{branch}.
            " which does not match other analysis/branch key:\t".$branch_codes->[0].' '.$branch);
        }  
      }
    }
  }
  else{#We only expect 1 branch numbers
    
    if(scalar(@$branch_codes) != 1){
      throw('Only 1 branch number is permitted if no branch config has been initialised, found '.
        scalar(@$branch_codes)); 
    }
    
    if($branch_codes->[0] !~ /[0-9]+/o){
      throw("Found analysis/branch key when no branch config has been initialised:\t".$branch_codes->[0].
            "\nPlease use a branch number or call init_branching_by_analysis or init_branching_by_config");
    }
    
    $branch = $branch_codes->[0]; 
  }

  return $branch;    
}



#$self->branch_job_group([['Branch2_analysis1', 'Branch2_analysis2'] [$branch2_job_id1, ...], 'RunIDR', [$job_id2, ...]]);
#will only ever have 1 semaphore

#This currently breaks the 'custom' branch, as we are validating all
#present in config. We could get around this, by default everything to the custom analysis
#if it is not present

sub branch_job_group{
  my ($self, $fan_branch_codes, $fan_jobs, $funnel_branch_code, $funnel_jobs) = @_;
  my $branch_config = $self->{branch_config}; #as this may be derived from analysis or method?
  my $fan_branch    = _get_branch_number($fan_branch_codes, $branch_config); 
  #this also asserts_ref for $fan_branch_codes
    
  if(! (check_ref($fan_jobs, 'ARRAY') &&
        scalar(@$fan_jobs) > 0)){
    throw('Must have at least 1 job in the job group');        
  }
  
  my $job_group = [$fan_branch, $fan_jobs]; 
   
  if(($funnel_branch_code && ! $funnel_jobs) ||
     ($funnel_jobs && ! $funnel_branch_code)){
    throw('Must have both a funnel branch code and funnel job(s) to define a funnel in a job group');    
  }
  elsif($funnel_branch_code){ #implicit that $funnel_jobs are also true
    #check $funnel_branch_code is not ref?
     
    if(! (check_ref($funnel_jobs, 'ARRAY') &&
          scalar(@$funnel_jobs) > 0)){
      throw('Must have at least 1 funnel job when defining a funnel in a job group'); 
    }
   
    my $funnel_branch = _get_branch_number([$funnel_branch_code], $branch_config); 
   
    if($branch_config){#check the funnel is valid wrt fan analyses
      
      #This will only work for branch config derived from the DataflowRules
      #not defined in the config!
      
      throw('How are we going to differentiate here?');
      
      if($branch_config->{$fan_branch_codes->[0]}{funnel} ne $funnel_branch_code){
        throw($funnel_branch_code.' is not a valid funnel analysis for '.$fan_branch_codes->[0].
          'Please check you dataflow configuration');  
      }    
    }
    
    push @$job_group, ($funnel_branch, $funnel_jobs);
  }
  
  $self->{job_groups} ||= [];
  push @{$self->{job_groups}}, $job_group;
  
  return;
}

#Woudl have been branched like this
#$self->branch_job_group([['Branch2_analysis1', 'Branch2_analysis2'], [$branch2_job_id1, ...], 'RunIDR', [$job_id2, ...]]);
#But should now look like this
##$self->branch_job_group([2, [$branch2_job_id1, ...], 3, [$branch3_job_id2, ...]]);

sub dataflow_job_groups{
  my $self       = shift;
  my $job_groups = $self->{job_groups} || []; #allow no job groups
  
  foreach my $job_group(@$job_groups){
     
    while(@$job_group){
      my $branch   = shift @$job_group;
      my $id_array = shift @$job_group;
      
      #let dataflow_output_id do the validation
      $self->dataflow_output_id($id_array, $branch);
    }
  } 

  return;
}

#we can't do branch_key_method dynamically from here
#as there is no generic way to support a method which
#get updated iteratively
#branch_key_method is validated in init_branch_config
#but has to be called prior top branch_ouput_id


#Haven't we given up on using a branch key method?
#check and remove!
#branch_key_param/method are to support
#generic analyses which may want to branch using a different param
#under different circumstances




=pod

sub branch_output_id{
  my ($self, $output_id, $branch_key, $branch) = @_;
  
  #These will fail ungracefully if we have not called init_branch_config
  #todo should we use AUTOLOAD to catch this and give a useful error about set/get_param_method not being called first?
  #can this be done just in Base?

  if(! (defined $branch ||
        defined $branch_key) ){
    throw('No valid branch_key_method or branch argument defined');  
  }
  
  #my $branch_output_ids = $self->branch_output_ids;
  
  #This is now resilient to a lack of branch_key config
  #but is this what we want?
  #do we want to force the use of init_branch_key_config, so
  #that we reduce the risk of missing it?
  
  
  #we probably don't want the branches hardcoding in the runnable as the dataflow
  #is really defined in the config. So let's keep the config there
  #and change the branch numbers to branch descriptors
  #This may seem like overkill, but it keeps analysis specific config
  #out of the runnable, and centralises the branching logic
  #hence reducing complexity in the Runnables themselves
  
   
  
  my $branch_output_ids = $self->get_param_method('branch_output_ids', 'silent', {});
  
  if($branch_key){
    my $branch_config     = $self->branch_config;
    $branch = (exists $branch_config->{$branch_key}) ? 
      $branch_config->{$branch_key} : 100;
    $self->helper->debug(2, "Branch key=$branch_key");   
  }
  
  if(! exists $branch_output_ids->{$branch}){
    $branch_output_ids->{$branch} = [];
  }
 
  $self->helper->debug(1, "Branching output_id to branch $branch:", $output_id);
  
  push @{$branch_output_ids->{$branch}}, $output_id;
  return;
}

### BRANCH_OTUPUT_ID/DATAFLOW_OUTPUT_ID ISSUES ###

### INITIAL REQUIREMENT #######
#
# 1 To allow defered data flow in the write_output method such that
#   we don't get unwanted jobs created if we specify -no_write
#   This was to support running IdentifySetInputs such that we can see what
#   would be kicked off without creating the jobs, in case our initial
#   set of parameters was incorrect. Alternative here would be to run a stand alone
#   job, but that doesn't have access to the pipeline config
#   THIS IS ESSENTIAL!
#   (Although we could just dataflow dependant on no_write?)

### ISSUES #######
#
# 1 We are adding config complexity, basically reversing the dataflow spec, 
#   omitting the semaphore details. This info is already available via the 
#   DataflowRuleAdaptor
#
# 2 Due to the lack of proper semaphore handling, the current methods are unsafe.
#   They will work fine if the config does not contain semaphores, but it currently
#   does check or handle grouping of fan/funnel jobs
#
# 3 Current implementation assumes highest branch number might be a funnel and so flows 
#   this last. This only works there is only ever 1 funnel job for the whole analysis and 
#   that funnel job happens to be the highest branch number.
#   This is non optimal in many way, not least because it runs counter to the idea of being able to 
#   dynamically expand the dataflow branches (dependant on adding analysis e.g. mutliple peak caller branches)
#   In this case, it is better to have the generic (in that they don't case which branch they are fed from,
#   e.g. can take input from RunCCAT or RunSWEmbl) funnels jobs first. Such that you can expand
#   the dataflow config with new peak caller branches, without the risk of having to change the branch
#   of the funnel analyses in the config or the runnable. e.g.
#
#      'A->2' => [ 'RunIDR' ],   #<-- this will never change
#      '3->A' => ['Preprocess_BWA_replicate'],
#      '4'    => ['Preprocess_BWA_merged'],     
#      '5'    => ['Preprocess_BWA_control'],  
#       #we can add new aligner preprocess jobs here
#      '6->A' => ['Preprocess_OTHERALIGNER_replicate'],
#      '7'    => ['Preprocess_OTHERALIGNER_merged'],     
#      '8'    => ['Preprocess_OTHERALIGNER_control'],  
#
# Note: it would be more safe to block these in sets of 10, such that we can increase the
# analysis flowed to for each analysis e.g.
#      '21->A' => ['Preprocess_BWA_replicate'],
#      '22'    => ['Preprocess_BWA_merged'],     
#      '23'    => ['Preprocess_BWA_control'],  
#      '24'    => ['SomeNew_BWA_analysis_here']
#
# Although this is merely for human readability, only having funnels first really matters
# in terms of maintainability. All the rest should be picked up and used correctly


### SOLUTIONS #####
#
# 1 Maintain branch_config and change branch_output_id to cache_job_group.
#   This will change from a simply hash cache to a 2d array of hashes.
#   e.g. $self->cache_job_group([{'Preprocess_BWA_replicate'=> $job_id1}, {'RunIDR' => $job_id2} ]);
#   The arrays are necessary to maintain the order of the dataflow for funnel jobs.
#   cache_job_group would push this array onto an array of job_groups 
#   Runnable will have to know which jobs need and associated funnel job.
#   This will not catch problems where job group do not have the required funnel
#   Nor will it catch job groups with mismatched funnel analyses
#   How are we going to support multiple jobs? 
#   This should look more like this:
#   e.g. $self->branch_job_group([{'Preprocess_BWA_replicate'=> [$job_id1, ...]}, {'RunIDR' => [$job_id2, ...]} ]);
#   what about sorting multiple imputs?
#   That's fine just define additional analyses as above, these can be validated vs the branch config


#
# 2 Drop branch_config and generate this via DataflowRuleAdaptor in init_branch config.
#   This will function as above, apart from that is iwll be posssible
#   to validate the funnel jobs are present and correct for a given job group
#   A cache should be built of valid funnel jobs for each branch name(fan analysis).
#   This may break if the Dataflow API/mechanism is changed

### CAVEATS #####
#
# 1 All this work on the assumption that each branch flow to only 1 analysis
#   i.e. the analysis we can build from the branck key value(aligner/peak caller name).
#   With solution 2, it is possible to check this and throw if >1 is configured or allow
#   if a allow_multiple_branch_analyses flag is passed.
#
# 2 All of these solutions require standardised branch name format based on the branch key
#
# 3 Funnel handling will always need to be hardcoded in the runnable, as it is impossible
#   to tell from the dataflow config should funnel after 1 branch input or all branch inputs
#   This is handled by the job grouping.
#
# 4 Due to caveat 3 funnel handling will be unsafe or impossible in runnables where the branch 
#   config and data flow can change between instances of that runnable as a different analysis
#   e.g. IdentifySetInputs
#   These analyses normally use a branch key method to define the branch name, and have no
#   inherant knowledge of dataflow or sempahores. i.e. they might be able to 
#   get the semaphore info for a given branch, but they might not know what to send down it
#   This would also be slightly problematic in building the job groups, conditionally
#   including the funnel job id.
#
# 5 This is really polluting the runnable with config info/requirements. Ideally we want them
#   to be entirely separated. However, this is the only way of doing the dynamic branching.
#   Hence moving this functionality as far away from the runnable as possible is preferable
#   In Base as present, but could be moved deeper into hive code with wrappers

### QUESTIONS
#
# 1 Where are semaphores defined in Dataflow rule? funnel_dataflow_rule_id
# 2 How does DataflowRule handle multiple target analyses? Are these separate rules?
# 3 Would need to inject DataflowRuleAdaptor::fetch_all_by_analysis_id
#   as it currently assumes you pass a branch code(or reserved name  MAIN, ANYFAILURE, RUNLIMIT, MEMLIMIT)
# 4 would also be good to inject get_analysis_dataflow_config which would return the config we require



sub dataflow_branch_output_ids {
  my $self = shift;
  my $dataflow_output_ids = $self->branch_output_ids;#set by init_branch_config and will always be defined
  
  if(keys %{$dataflow_output_ids}){    
    #We need to flow the job in the order they are specified in the config
    #i.e. we can't flow the accumulator(branch 1) before we have flown the fan jobs
    #hence descending sort
        
    foreach my $branch(sort {$b <=> $a} (keys %$dataflow_output_ids)){
     
      $self->helper->debug(1, 'Dataflowing '.scalar(@{$dataflow_output_ids->{$branch}}).
        " job(s) on branch $branch");    
      $self->dataflow_output_id($dataflow_output_ids->{$branch}, $branch);
    }
  }  
  
  return;
}

=cut 

#There is effectively no difference between these apart from
#dataflow_params are meant to be specified at different levels in the config
#dataflow_params are to be specified in analysis parameters for flowing only 
#onto the next job
#batch_params are to be defined in the pipeline_wide_params and are to be flow
#for all appropriate analysis for that batch. It might not be appropriate
#to flow these for some fan jobs, but this can be handled in the analysis on a case
#by case basis

#dataflow_params used to be defined in the input_id and be completely flexibly
#e.g.{dataflow_params=>{recover=>1}}
#It is still possible to over-ride the harcoded defaults
#but it must be done as follows:
#{dataflow_params=> ['recover'], recover=>1}
#This also makes the params available to the intial analysis
#rather than having them nested in the dataflow_params hash


#These are only dataflowed to the next analysis

sub dataflow_params {
  my $self         = shift;
  my $optional     = shift;
  my $param_values = {};
  
  if(! $self->param_silent('dataflow_param_names')){
    
    if(! $optional){
      throw('dataflow_param_names are required but not defined.'.
        ' Please define or set the optional flag');  
    }
  }
  else{
    $param_values = $self->_dataflow_params_by_list('dataflow_param_names') 
  }
  
  return $param_values;
}


sub batch_params {
  return shift->_dataflow_params_by_list('batch_param_names');
}


sub _dataflow_params_by_list {
  my ($self, $method_list_name) = @_;  
  
  #Can redefine params inline by simply redefining the param in question
  #This will however change the batch_param for any
  #subsequent jobs.
  
  #This should only ever be called once by the Runnable in question
  #so let's not cache
  #this will also allow dynamic updating of batch params
  #should one of them change between calls of batch_params method
    
  my %param_values = ();#Define as we will deref in the caller
  
  $self->helper->debug(3, "Getting $method_list_name ", $self->param_required($method_list_name)); 
   
  foreach my $param_name(@{$self->param_required($method_list_name)}){
  
    my $param = $self->param_silent($param_name);
    
    if( blessed($param) ){
      
      if(! exists $object_dataflow_methods{ref($param)}){
        throw("Cannot dataflow $method_list_name $param_name as it is an object without a defined dataflow method:\t".
          $param."\nPlease add to \%object_dataflow_methods");  
      }
      else{
        my $data_flow_method = $object_dataflow_methods{ref($param)};
        $param = $param->$data_flow_method;
      }
    }
    
    if(defined $param){
      $param_values{$param_name} = $param;
    }
  }
  
  return \%param_values;
}


#Could also extend the batch_param_names for the rest of this branch of the pipeline

#Changed from check_link_analysis_can_run, so is now generic.
#This would also provide a simply mechanism for stopping a pipeline at a given point
#without borking the potential ongoing dataflow
#This could be useful for various reasons.
#although doing this in a batch manner may be problematic, as resetting this variable to 1
#would require fiddling with the individual job IDs

#To fail or not to fail?
#1 Fail. Set transient_error to 0 and die with a useful error message.
#  This will mask true failures for other reasons. But the job will still be on
#  there for hive to pick up, once we have updated the reelvant can_run_AnalysisLogicName param
#  and reset the job. Reset, will be to remove the transient error to 1 (and the retry count?)
#  We will be able to manage this somewhat by having a function to discriminate between true failures
#  and those triggered as above. This will involve checking the transient_error status and pattern matching
#  the log message.
#2 Succeed but do nothing. This will not hide true failures and will also reflect the true status of the
#  job as defined by our criterion i.e. we want it to succeed in doing nothing.
#  This however makes if very hard to identify outstanding jobs which have been halted. This means
#  that when continued dataflow is required(after topup or whatever), then it is almost impossible
#  to identify those job which were halted. Hence these jobs would need to be reseeded. This is non-trivial
#  as we may not know the original seed spec, and it is likely to create duplicate job ids, which would 
#  then do nothing if the can_run_AnalsisLogicName param is not in the job id (i.e. it is pipeline_wide).
 
#param_require should/will set transient_error to 0 before dying if not specified
#as opposed to specified as undef. Hence reducing retry churn.

#Fail or no fail flag and return boolean?

sub check_analysis_can_run{
  my $self      = shift; 
  my $check_int = $self->param_silent('check_analysis_can_run');
  my $can_run   = 1;
  
  if(defined $check_int){ #Might be undef, in which case we can run
    my $lname     = $self->analysis->logic_name;
    my $run_param = 'can_'.$lname; 
    $can_run      = $self->param_required($run_param);
   

    if(! $can_run){ #0 or specified as undef
      $self->input_job->transient_error(0); #So we don't retry  
      die("$lname has aborted as $run_param == 0 || undef"); 
      #Parse this in the env
      #Omit from GetFailedJobs 
      #and add GetAbortedJobs?
    }  
  }
    
  return;
}



#todo Move these to (and create) BaseSequenceAnalysis.pm?
#We really need to specify a force_male flag if gender is undefined

#species and assembly are set in _set_outdb

sub sam_ref_fai {
  my $self   = shift;
  my $gender = shift; 
  
  if(! defined $self->param_silent('sam_ref_fai')){

    if(! defined $gender){
      $gender = $self->param_silent('gender') || $self->param_silent('default_gender');
      
      if(! defined $gender){
        $self->throw_no_retry('No gender argument or param defined and no default_gender '.
        'specific in the config');
      }
    }
    
    my $file_name = $self->species.'_'.$gender.'_'.$self->assembly.'_unmasked.fasta.fai';
    my $sam_ref_fai = validate_path([$self->data_root_dir,
                                    'sam_header',
                                    $self->species,
                                    $file_name]);
     
    $self->param('sam_ref_fai', $sam_ref_fai);
  }
  
  return $self->param('sam_ref_fai');
}

#merge these two subs

sub sam_header{
  my $self   = shift;
  my $gender = shift; 
  
  if(! defined $self->param_silent('sam_header')){
  
    if(! defined $gender){
      $gender = $self->param_silent('gender') || $self->param_silent('default_gender');
      
      if(! defined $gender){
        $self->throw_no_retry('No gender argument or param defined and no default_gender '.
        'specific in the config');
      }
    }
    
    my $sam_header = join('/', ($self->param_required('data_root_dir'),
                                 'sam_header',
                                 $self->species,
                                 $self->species.'_'.$gender.'_'.$self->assembly.'_unmasked.header.sam'));
    $self->set_param_method('sam_header', $sam_header);
  }
  
  return $self->param('sam_header'); 
}



sub get_alignment_file_prefix_by_ResultSet{
  my ($self, $rset, $control) = @_;  
    
  $self->out_db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::ResultSet', $rset);
  my $path = $self->alignment_dir($rset, undef, $control);
  my @rep_numbers;
  
  foreach my $isset(@{$rset->get_support}){
    
    if( ($control && $isset->is_control)  ||
        ((! $control) && 
         (! $isset->is_control)) ){ 
       push @rep_numbers, $isset->replicate;
    }
  }
  
  #Uses get_set_prefix_from_Set rather than set name, as IDR rep sets
  #will already contain the rep number, and merged sets will not  
  return $path.'/'.$self->get_set_prefix_from_Set($rset, $control).'_'.
    $rset->analysis->logic_name.'_'.join('_', sort(@rep_numbers)); 
}

#This will have validated that input_set are all part of the same experiment
#(and they have the same alignment logic_name, when this is implemented) ???

sub get_alignment_files_by_ResultSet_formats {
  my ($self, $rset, $formats, $control, $all_formats, $filter_format) = @_;
  assert_ref($formats, 'ARRAY');  
  my $file_type = ($control) ? 'control_file' : 'alignment_file';
  my ($path, $align_files);

  #Currently don't support > 1 result_file, but that is caught in the Importer somewhere
  #Move this logic to the Importer and make input_files/data_dir optional
  #Allow file over-ride here
  #usage of get_param_method here will mean we can't call this again
  #for the same file type
  #hence, this will likely fail if the worker runs more than 1 job
  
  
  if(! $rset->has_status('ALIGNED')){
    throw("Cannot get alignment files for ResultSet which does not have ALIGNED status:\t".
      $rset->name);
  }
   
  if($self->get_param_method($file_type, 'silent')){ #Allow over-ride from config/input_id   
    #Need to test in here that it matches one of the formats
    #Need to add support for converting this non-standard path to other formats
    throw("$file_type config over-ride is not yet implemented");
    #$align_files = { => validate_path($self->param($file_type)) };
  }
  else{ # Get default file
    my $params = {ref_fai            => $self->sam_ref_fai,  #Just in case we need to convert
                  filter_from_format => $filter_format,
                  all_formats        => $all_formats};  
    
    $path = $self->get_alignment_file_prefix_by_ResultSet($rset, $control).
                '.samse'; #Currently hardcoded for bam origin!
  
    $path .= '.unfiltered' if $filter_format;
     

    $self->helper->debug(1, "Getting $file_type (filter_from_format=$filter_format):\n\t".$path);
    $align_files = get_files_by_formats($path, $formats, $params);
  }  
 
  #throw here and handle optional control file in caller. This should be done with 
  #a no_control/skip_control flag or similar  
  return $align_files || throw("Failed to find $file_type (@$formats) for:\t$path");  
}


sub validate_package_from_path{
  my $self     = shift;
  my $pkg_path = shift;  
  
  eval { require $pkg_path; };

  if ($@) {
    throw( "Failed to require:\t$pkg_path" );
  }
  
  #This might not always be correct if the file contains >1 package
  return path_to_namespace($pkg_path);
}


#config. broad_peak_feature_types is now required, even if it is just and empty array
#as this will raise awareness that the pipeline is trying to perform some idr
#sensitive functions
  
#This ! grep is fine, although it returns an empty string instead of 0
#as oppose to 1, when nothing is returned from grep

sub is_idr_feature_type{
  my $self  = shift;
  my $ftype = shift;
  throw('Must pass a FeatureType name') if ! defined $ftype; 
  
  return $self->param_silent('no_idr') ? 0 : 
          ! grep(/^${ftype}$/, @{$self->param_required('broad_peak_feature_types')});
}



#Move this to EFGUtils

=head2 C<CvGV_name_or_bust> I<coderef>

Calls L<Devel::Peek> to try to find the glob the ref lives in; returns
C<undef> if L<Devel::Peek> can't be loaded, or if C<Devel::Peek::CvGV> can't
find a glob for this ref.

Returns C<< I<package>::I<glob name> >> if the code ref is found in a glob.

=cut

# From the perl debugger
# Could alternatively use the more simple Sub::Identify implementation 
# but this is not a core module

sub CvGV_name_or_bust {
    my $in = shift;
    return unless ref $in;
    $in = \&$in;            # Hard reference...
    #eval { require Devel::Peek; 1 } or return;Now used above
    my $gv = Devel::Peek::CvGV($in) or return;
    *$gv{PACKAGE} . '::' . *$gv{NAME};
} ## end sub CvGV_name_or_bust



# Put this in Process.pm
  
=head2 throw_no_retry

  Arg [1]    : string $msg
  Arg [2]    : (optional) int $level
               override the default level of exception throwing
  Example    : $self->throw_no_retry('We have a non transient error, this job will not be retried');
  Description: Sets transient_error to 0 before throwing, such that this job will not be retried
  Returntype : none
  Exceptions : throws
  Caller     : general

=cut


sub throw_no_retry {
  my $self = shift;
  $self->input_job->transient_error( 0 );
  throw(@_);

} ## end throw_no_retry  

1;
