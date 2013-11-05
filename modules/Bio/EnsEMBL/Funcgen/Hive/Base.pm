=pod 

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
                                               

use base ('Bio::EnsEMBL::Hive::Process');

#todo what won't be using the tracking DB? i.e. Isn't everything a BaseDB?
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

sub alignment_dir {
  my ($self, $input_set, $create) = @_;
   
  if(defined $input_set){ 
    my $exp_id = $input_set->get_Experiment->dbID;
    
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
                                 $self->get_study_name_from_InputSet($input_set)],
                                $create);    
  }
  return $self->param('alignment_dir');
}


sub get_study_name_from_InputSet {
  my $self = shift;
  my $iset = shift;
      
  return $self->parse_study_name($iset->get_Experiment->name,
                                 $iset->cell_type->name,
                                 $iset->feature_type->name);
}


#Just strip off first two words of experiment name for now
#until we model study/experiment properly wrt having multiple input_sets associated
#This should still work if the experiment names have ctype and ftype removed

sub parse_study_name_name{
  my ($self, $exp, $ctype, $ftype) = @_;   
  (my $study_name = $exp) =~ s/${ctype}_${ftype}_(.*)/$1/;
  return $study_name;
}


sub get_output_work_dir_methods{
  my ($self, $default_odir, $no_work_dir) = @_;
   
  my $out_dir = $self->validate_dir_param('output_dir', 1, $default_odir); 
   
   
  my $work_dir;    
             
  if(! $no_work_dir){                 
    my $dr_dir             = $self->data_root_dir;
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

#we dei=finitely need branch_key_param|method
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

sub init_branch_config {
  my ($self, $optional, $validate_bkey_method) = @_;  
  my $branch_config = $self->get_param_method('branch_config', 'silent');
 
  if(! $optional){
    assert_ref($branch_config, 'HASH', 'branch_config');  
  }
  
  if( (defined $branch_config) &&
      $validate_bkey_method){
    my $bkey_method = $self->get_param_method('branch_key_method', 'required');
    #don't ever need this param again, so don't get_param_method
    
      
    #warn "branch_key_method is $bkey_method";  
      
    #This will use an existing method or create one from en exiting parameter
    
    if(! defined $self->get_param_method($bkey_method, 'required')){
      throw("Required branch_key_method $bkey_method is not defined defined");
    }
    
    #Now reset branch_key_method method to actual branch_key_method?
    #or just go down the old method as a variable root  
  }
  
  $branch_config ||= {};
  
  #100 is reserved branch for custom analyses
  my %branch_data_flow_ids = ();
  #don't set 1 and 100 now so we don't have
  #keys if we really don't flow anything
  map { $branch_data_flow_ids{$_} = [] } values %$branch_config;
  
  $self->helper->debug(1, "Initialised branch_data_flow_ids:", \%branch_data_flow_ids);
  
  $self->set_param_method('branch_output_ids', \%branch_data_flow_ids);
    
  return;
}

#we can't do branch_key_method dynamically from here
#as there is no generic way to support a method which
#get updated iteratively
#branch_key_method is validated in init_branch_config
#but has to be called prior top branch_ouput_id

sub branch_output_id{
  my ($self, $output_id, $branch_key, $branch) = @_;
  
  #These will fail ungracefully if we have not called init_branch_config
  #todo should we use AUTOLOAD to catch this and give a useful error about set/get_param_method not being called first?
  #can this be done just in Base?

  if(! (defined $branch ||
        defined $branch_key) ){
    throw('No valid branch_key_method or branch argument defined');  
  }
  
  my $branch_output_ids = $self->branch_output_ids;
  my $branch_config     = $self->branch_config;
  
  if($branch_key){
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


sub dataflow_branch_output_ids {
  my $self = $_[0];
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
  my ($self, $optional) = @_; 
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
  return $_[0]->_dataflow_params_by_list('batch_param_names');
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



#todo Move these to BaseSequenceAnalysis?

sub sam_ref_fai {
  my $self = $_[0]; 
  
  if(! defined $self->param_silent('sam_ref_fai')){
    my $file_name = $self->species."_"; 
    $file_name  .= $self->param_required('cell_type')->gender || 'male'; #why isn't this female?
    $file_name  .= '_'.$self->assembly.'_unmasked.fasta.fai';
    
    my $sam_ref_fai = validate_path([$self->data_root_dir,
                                    'sam_header',
                                    $self->species,
                                    $file_name]);
     
    $self->param('sam_ref_fai', $sam_ref_fai);
  }
  
  return $self->param('sam_ref_fai');
}



#Should really also pass filter_from_format
#as we might not always want to filter
#and we may pick up the original filter_from_format config
#so expicilty pass it rather than messing around with redefining and dataflowing
#This won't fix the PreprocessAlignments issue!
#As CollectionWriter is used again.


#alignment will use this to generate the output file

#This does not account for replicates?
#Replicate naming is handling via input set name for signal sets
#control files will always be merged and so don't need the TR suffix
#from the input set name

sub get_alignment_file_prefix_by_InputSet{
  my ($self, $iset, $control) = @_;  
  
  $self->out_db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::InputSet', $iset);
  my $path = $self->alignment_dir($iset);
  
  if($control){
    $path .= '/'.$iset->cell_type->name.'_'.
        $self->param_required('control_feature').'_'.
        $self->get_study_name_from_InputSet($iset);
  }
  else{
    $path .= '/'.$iset->name; 
  } 
  
  return $path; 
}

#This will have validated that input_set are all part of the same experiment
#(and they have the same alignment logic_name, when this is implemented) ???

sub get_alignment_files_by_InputSet_formats {
  my ($self, $iset, $formats, $control, $all_formats, $filter_format) = @_;
  assert_ref($formats, 'ARRAY');  
  my $file_type = ($control) ? 'control_file' : 'alignment_file';
  my ($path, $align_files);

  #Currently don't support > 1 result_file, but that is caught in the Importer somewhere
  #Move this logic to the Importer and make input_files/data_dir optional
  #Allow file over-ride here
  #usage of get_param_method here will mean we can't call this again
  #for the same file type
  #hence, this will likely fail if the worker runs more than 1 job
  
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
    
    $path = $self->get_alignment_file_prefix_by_InputSet($iset, $control).
                '.samse'; #Currently hardcoded for bam origin!
  
    if($filter_format){
      $path .= '.unfiltered';  
    }

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


  
  

1;
