=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( get_study_name_from_Set
                                               get_set_prefix_from_Set
                                               scalars_to_objects 
                                               validate_path
                                               run_system_cmd
                                               dump_data );
use Bio::EnsEMBL::Funcgen::Hive::Utils     qw( inject_DataflowRuleAdaptor_methods );
use Bio::EnsEMBL::Funcgen::Sequencing::SeqTools qw( get_files_by_formats );
use Bio::EnsEMBL::Utils::Scalar            qw( assert_ref check_ref );  
use Scalar::Util                           qw( blessed );

use base qw(Bio::EnsEMBL::Hive::Process);

my %debug_modes = (no_tidy   => 1,
                   no_tidy_1 => 1,
                   no_tidy_2 => 2,
                   no_tidy_3 => 3);

#Used in set_param_arrays for scalars_to_objects
#Could our this to make it dynamically mutable, i.e. redefine it from a sub class
#Need this here rather than EFGUtils so we can test
my %param_class_info = 
 (epigenomes          => ['Epigenome',         'fetch_by_name'],
  feature_types       => ['FeatureType',       'fetch_by_name'],
  analyses            => ['Analysis',          'fetch_by_logic_name'],
  experimental_groups => ['ExperimentalGroup', 'fetch_by_name'],
  experiments         => ['Experiment',        'fetch_by_name'],
  epigenome           => ['Epigenome',         'fetch_by_name'],
  feature_type        => ['FeatureType',       'fetch_by_name'],
  analysis            => ['Analysis',          'fetch_by_logic_name'],
  experimental_group  => ['ExperimentalGroup', 'fetch_by_name']       );


my %object_dataflow_methods = ('Bio::EnsEMBL::Analysis' => 'logic_name');

# my %valid_file_formats = 
#  (bed  => 'Bed',
#   sam  => 'SAM',
#   bam  => 'BAM' );

#Advantage of having separate fetch_input, run and write methods is to allow
#calling of super methods at appropriate point.
#Theoretically, this could all be done in the run method with calls to the super methods at the 
#appropriate points, but this would break the object model and make it unihneritable.

#This defines a set of parameters based on given parameters from the pipeline:
sub fetch_input {   # nothing to fetch... just the DB parameters...
  my $self = $_[0];
  $self->SUPER::fetch_input;
   
  if(my $debug_level = $self->debug){
    #This debug mode handling doesn't work yet, as the main hive code barfs 
    #if the value is not a number
    #let's get leo to look at this, as it is quite useful
  
    if($debug_level !~ /^[1-3]$/){
      
      if(! exists $debug_modes{$debug_level}){
        $self->throw("Not a valid -debug mode:\t$debug_level\nPlease specify one of:\t".
          join(' ', keys %debug_modes));  
      }
      
      (my $debug_mode = $debug_level) =~ s/_[0-3]$//o;
      $self->param($debug_mode, $debug_modes{$debug_level});
    }
  
    $main::_debug_level = $debug_level; #Pass this onto the Helper and any other modules in the process
  } 
  #Catch and validate generic params here
  #but manage mandatory aspect in more specific runnables
  
  #validate_dir_params uses process_params which uses get_param_and_method
  #Hence all these dir with now have a method available
  
  #root_input_dir is mandatory and workdir and hive_output_dir have defaults based on this
  #Mandatory
  my $data_root_dir = $self->param('data_root_dir');
  if (! -d $data_root_dir) {
    system("mkdir -p $data_root_dir");
  }
  
  $self->validate_dir_param('data_root_dir'); 
 
  $self->validate_dir_param('bin_dir',    undef, 1); #optional

  $self->validate_dir_param
   ('work_root_dir', 
    1,
#     $self->data_root_dir.'/output/'.$self->param_required('pipeline_name')
    undef
   );

  #output and work dir need setting based on dbname
  if($self->get_param_method('use_tracking_db')){
    $self->_set_out_db;      
  }

  #Sneakily tidy up any tmp files defined by the previous job
  #This is normally a funnel job, and the tmp files
  #were generated in the job which created the fan/semphore. 
  #The tmp files will have been used in the fan jobs . 
  #It is unsafe to delete them in the fan jobs in case they 
  #die (RUN_LIMIT) before the job is marked as DONE
  #We only ever want to get this param (not_set), as would always pass this explicitly in the 
  #relevant output_id and not generically to all branches via dataflow_params or batch_params
  my $garbage = $self->param_silent('garbage');


  #Should garbage collection and archiving be disabled if -no_write is set?

  if(defined $garbage){
    #allow scalars and arrayref
    
    if(ref($garbage)){
      assert_ref($garbage, 'ARRAY', 'garbage files');
      unlink(@$garbage);
    }
    else{
      unlink($garbage);  
    }    
  }
  
  my $to_archive = $self->param_silent('to_archive');
  $self->archive_files($to_archive) if defined $to_archive;
    
  return;
}

# Validate mandatory non DB params separately 
# as they are not mandatory for subclass BaseDB
# This should be called for any cubclass which is not BaseDB

sub validate_non_DB_inputs{
 my $self = $_[0];
 
 $self->get_param_method('species', 'required');
 $self->get_param_method('assembly', 'required'); 
 return;
}

sub hive_run_system_cmd {

  my $self    = shift;
  my $command = shift;
  my $no_exit = shift;
  my $verbose = shift;
  
  print "$command\n" if $verbose;
  system($command);
  return $self->_hive_handle_exit_status($command, $?, $!, $no_exit);
}

sub _hive_handle_exit_status {

  my $self         = shift;
  my $cmd          = shift;
  my $exit_status  = shift;
  my $errno        = shift;
  my $no_exit      = shift;
  
  my ($exit_code, $err_string);

  $exit_code = $exit_status >> 8;

  if ($exit_status == -1) {
    $err_string = "Failed to execute:\t$cmd\nError:\t$errno\n";
  }
  elsif ($exit_status & 127) {
    $err_string = sprintf("Child process died with signal %d, %s coredump\n$cmd\nError:\t$errno\n",
    ($exit_status & 127),
    ($exit_status & 128) ? 'with' : 'without');
  } elsif ($exit_status != 0) {
    $err_string = sprintf("Child process exited with value %d:", $exit_code)."\t$cmd\nError:\t$errno\n";
  }

  if(defined $err_string) {
    if(! $no_exit) {
      $self->throw($err_string);
    } else {
      warn $err_string;
    }
  }
  return $exit_code;
}

sub alignment_root_dir {
 my $self = shift;
 
 if(! $self->param_silent('alignment_root_dir')) { 
   $self->set_dir_param_method(
      'alignment_root_dir', [
	$self->bam_output_dir,
      ], 1
    );
 }
 return $self->param('alignment_root_dir'); 
}

sub alignment_dir {

  my ($self, $rset, $create, $control) = @_;
  
  if(defined $rset) {

    $self->set_dir_param_method(
      'alignment_dir',
      [
	$self->alignment_root_dir,
 	get_study_name_from_Set($rset, $control)
      ],
      $create
    );
  }
  
  return $self->param('alignment_dir');
}

use File::Spec;

=head1 regulation_directory

  The directory in which the regulation files go. This is shared with the 
  other teams, so best to put it in its own directory.
  
  Deactivated for now, use own directory after the development break when
  we move from result_set to alignment.

=cut
sub regulation_directory {
  return undef;
  return 'funcgen';
}

=head1 version_directory

  The Ensembl release version when the files were generated. Padded with a 
  zero so it looks nice when we hit release 100.

=cut
sub version_directory {
  my $self = shift;
  return $self->param('ensembl_release_version');
}

sub db_output_dir {
  my $self = shift;
  return File::Spec->catfile(
    $self->data_root_dir,
    $self->species,
    $self->assembly,
    &regulation_directory,
  );
}

=head1 default_directory_by_table_and_file_type

  Returns the names of directories in which the files should be stored by 
  default. This is:
  
  <root directory> / <the directory forregulation specific files> / <the name of the table that represents it in the database> / < the release version > / ersa_signal / <the file type>

=cut
sub default_directory_by_table_and_file_type {
  my $self      = shift;
  my $table     = shift;
  my $file_type = shift;
  
  return File::Spec->catfile(
    $self->db_output_dir, $table, $self->version_directory, 'ersa_signal', $file_type
  );
}

sub peaks_output_dir  {  return shift->default_directory_by_table_and_file_type('annotated_feature', 'peaks');  }
sub bam_output_dir    {  return shift->default_directory_by_table_and_file_type('result_set',        'bam');    }
sub bigwig_output_dir {  return shift->default_directory_by_table_and_file_type('result_set',        'bigwig'); }

=head1 quality_check_output_dir

  The root directory for the output of all quality checks.

=cut
sub quality_check_output_dir {
  my $self      = shift;
  return File::Spec->catfile(
    $self->db_output_dir, 'quality_checks', $self->version_directory
  )
}

sub fastqc_output_dir                       { File::Spec->catfile( shift->quality_check_output_dir, 'fastqc')                       }
sub flagstats_output_dir                    { File::Spec->catfile( shift->quality_check_output_dir, 'flagstats')                    }
sub phantom_peaks_output_dir                { File::Spec->catfile( shift->quality_check_output_dir, 'phantom_peaks')                }
sub proportion_of_reads_in_peaks_output_dir { File::Spec->catfile( shift->quality_check_output_dir, 'proportion_of_reads_in_peaks') }
sub chance_output_dir                       { File::Spec->catfile( shift->quality_check_output_dir, 'chance')                       }

sub temporary_directory_root {
  return File::Spec->catfile(
    shift->work_root_dir, 'temp'
  )
}

sub get_output_work_dir_methods {

  my $self         = shift;
  my $default_odir = shift;
  my $no_work_dir  = shift;
  my $out_dir      = $self->validate_dir_param('output_dir', 1, $default_odir);  # Create flag
  my $work_dir;    
             
  if(! $no_work_dir){
    my $dr_dir = $self->data_root_dir;
    
    if($out_dir !~ /$dr_dir/){
      $self->throw('Cannot set a work_dir from an output dir('.$out_dir.
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

  $default_odir ||='';#to avoid undef warning in debug
  $self->helper->debug(1, "output_dir $out_dir\n\twork_dir $work_dir\n\tdefault output_dir$default_odir");
  return ($out_dir, $work_dir);
}

sub no_write     { return $_[0]->param_silent('no_write');  } 

sub _set_out_db {
  my $self = $_[0];
  
  #This will either be the params hash or a DBADaptor
  my $db = $self->get_param_method('out_db', 'required');
 
  if(! ref($db)){ #Has been set to scalar in config?
     $self->throw("out_db config is not a Hashref of DBAdaptor::new parameters:\t$db");
  }
  #Always required as we don't want to default to ensembldb for a pipeline!
  my $dnadb_params = $self->param_required('dnadb');  
      
  #Create TrackingAdaptor here, as we can't get_TrackingAdaptor later
  my $adaptor_class = ($self->use_tracking_db) ? 
    'Bio::EnsEMBL::Funcgen::DBSQL::TrackingAdaptor' 
    :
    'Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor'
    ;
    
  eval {
    $db = $adaptor_class->new(%{ $db }, %{ $dnadb_params });
  };
  if($@) {
    $self->throw("Error creating the $adaptor_class\n$@");  
  }
  if(! $db->isa('Bio::EnsEMBL::DBSQL::BaseAdaptor')) {
    $self->throw("The out_db param is set to an unexpected reference:\t" . (ref $db) . "\n"
      . Dumper($db)
      . "Please check the out_db config, which should be a Hashref of DBAdaptor parameters");
  } 

  #Reset the $db to an actual DBAdaptor (rather than a Tracking/BaseAdaptor)
  if ($self->use_tracking_db){
    $self->set_param_method('tracking_adaptor', $db);
    $db = $db->db; 
  }
    
  # Test connections
  $db->dbc->db_handle->ping();
  $db->dnadb->dbc->db_handle->ping();

  $db->dbc->disconnect_when_inactive(1);
  $db->dnadb->dbc->disconnect_when_inactive(1);
  
  # Skip, if module is run as standalone job.
  if ($self->dbc) {
    $self->dbc->disconnect_when_inactive(1);
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
      $self->throw("Failed to identify default assembly from the DB");
    }
    
    $self->assembly($assembly);
  }
  else{ #validate it exists in the dnadb
      
    my $cs = $cs_adaptor->fetch_all_by_version($assembly)->[0];
  
    if(! defined $cs){
      $self->throw("The -assembly $assembly does not exist in the database:\t".
	$db->dnadb->dbc->dbname); 
    }
    #todo make sure all Slice/CoordSystem calls use this assembly explicitly.
  }
  
  # VALIDATE/SET Species
  my $param_species = $self->get_param_method('species', 'silent');
  my $db_species = $db->species;
  
  if( ($param_species && $db_species) &&
    (lc($param_species) ne lc($db_species)) ){
      $self->throw('Mismatch between the DB species (meta species.production_name) and'.
      " -species parameter:\t".$db_species."\t".$param_species);      
  }
  elsif(! defined $param_species){
  
    if(! defined $db_species){
      $self->throw('Must either set -species or use a DB where the correct species '.
	'name is set in the meta table as species.production name');  
    }
  
    # Have to be mindful that we can't access the adaptors from the 
    #Registry using the -species parameter. Hence it is unsafe to use the Registry!
    #Always get the adaptors directly from the DB.
    $self->species($db_species);
  }


  
  return $self->out_db($db);
}



#Expose param_silent for handling optional params

sub param_silent {
  my $self = shift;
  return $self->input_job->_param_silent(@_); 
}

sub process_params {
  my ($self, $param_names, $optional, $as_array) = @_;
  
  #handle single scalar param name
  my $ref_type = ref($param_names);
  if($ref_type){ #can be ''
    
    if($ref_type ne 'ARRAY'){
      $self->throw("Param names argument cannot be a $ref_type.\n".
        'Only a single scalar or an Arrayref of param names is permitted.');  
    }
    
    if(! scalar @$param_names){
      $self->throw("Cannot pass an empty Arrayref for param names argument\n".
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
        $self->throw('process_params can only handle scalar (comma separated), '.
          'Arrayref of scalars or Hashref values');  
      } 
      
      #warn "before param_class_info test with $param_name";
             
      if(exists $param_class_info{$param_name}){
 
        if($ref_type eq 'HASH'){#can only be hash
          $self->throw("Need to implement $param_name fetch_all for process_params");
        }
        else{ #We have an Arrayref of scalars
          $values = scalars_to_objects($self->out_db, @{$param_class_info{$param_name}},
                                       $values);
        }
        
        $self->helper->debug(2, "$param_name is now ", $values);
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
  
  
  $self->helper->debug(1, 'Processed param names('.join(' ', @$param_names).
    ") are:\t".dump_data(\%all_params));
  
  return \%all_params; 
}

sub validate_dir_param {
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
  
  if($path){
  
    if($path =~ /.+\/$/o || $path =~ /\/\//o){
      $path =~ s/\/$//o;
      $path =~ s/\/\//\//o;
      $self->$dir_name($path);  
    }
   validate_path($path, $create, 1, $dir_name);
  }
  return $path;
}

sub set_dir_param_method {
  my ($self, $dir_param_name, $path, $create) = @_;      
  
  if(! (defined $dir_param_name &&
        defined $path) ){
    $self->throw('You must provide directory parameter name and path arguments to validate');
  }

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
      $self->throw("$req is not a valid method type for set_param_and method, can only be 'required'");  
    }
    elsif(! defined $param_value){
      # why not leave this to _param_and_method?
      $self->throw($param_name.' value is required but not defined');  
    }
  }
  
  return $self->_param_and_method($param_name, $param_value, $req);    
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
        $self->throw("Caught attempt to redefine non-Funcgen method:\t$package_method\n".
          "Please use an alternative method name, or use the above method directly if appropriate");
      }
    }
  }

  my $param_method = 'param';
  
  if(! defined $param_value) {
  
    if(defined $req_or_silent) {
    
      if( ($req_or_silent ne 'silent' ) &&
          ($req_or_silent ne 'required') ){
        $self->throw("Param method type can only be 'silent' or 'required', not:\t$req_or_silent");    
      }
    
      $param_method .= '_'.$req_or_silent;
    }
  }
  
  if(! defined $cref){
    #Inject method, hopefully not Pragha Khan style!
    no strict 'refs';
  
    *{ref($self)."::${param_name}"} = sub { 
      my $self  = shift;
      my $param = shift;
    
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
  
sub _get_branch_number{
  my $self          = shift;
  my $branch_codes  = shift;
  my $branch_config = shift;  
  assert_ref($branch_codes, 'ARRAY', 'Branch code array');
  my $branch;
 
  if(defined $branch_config){
    
    foreach my $bcode(@$branch_codes){
      
      $self->throw_no_retry('Branch code is undefined') if ! defined $bcode;
      #This will allow null string and 0
      
      if(! exists $branch_config->{$bcode}){
        # numbered or named branches might not actually be wired
        # dependant on loaded config
        warn "Could not find branch config for branch:\t".
            $bcode."\n";        
      }
      else{
        $branch ||= $branch_config->{$bcode}{branch};
        
        if($branch_config->{$bcode}{branch} != $branch){
          $self->throw($bcode.' fan analysis/branch key has branch '.$branch_config->{$bcode}{branch}.
            " which does not match other analysis/branch key:\t".$branch_codes->[0].' '.$branch);
        }  
      }
    }
  }
  else{#We only expect 1 branch number
    
    if(scalar(@$branch_codes) != 1){
      $self->throw('Only 1 branch number is permitted if no branch config has been initialised, found '.
        scalar(@$branch_codes)); 
    }
    
    if($branch_codes->[0] !~ /[0-9]+/o){
      $self->throw("Found analysis/branch key when no branch config has been initialised:\t".$branch_codes->[0].
            "\nPlease use a branch number or call init_branching_by_analysis or init_branching_by_config");
    }
    
    $branch = $branch_codes->[0]; 
  }

  return $branch;    
}

sub branch_job_group{
  my ($self, $fan_branch_codes, $fan_jobs, $funnel_branch_code, $funnel_jobs) = @_;
  my $branch_config = $self->{branch_config}; #as this may be derived from analysis or method?
 
  $fan_branch_codes = [$fan_branch_codes] if ! ref($fan_branch_codes);
  my $fan_branch    = $self->_get_branch_number($fan_branch_codes, $branch_config); 
  #this also asserts_ref for $fan_branch_codes
  
  if (!defined $fan_branch) {
    use Carp;
    use Data::Dumper;
    confess(
      "Can't find fan branch for (" 
      . Dumper( [ $fan_branch_codes, $branch_config ] )
      . ")");
  }
    
  $self->helper->debug(1, "Branching ".scalar(@$fan_jobs).
    " jobs to branch(es) $fan_branch(codes=".join(' ', @$fan_branch_codes).')');
    
  if(! (check_ref($fan_jobs, 'ARRAY') &&
        scalar(@$fan_jobs) > 0)){
    $self->throw("Must have at least 1 job in the job group for fan/branch $fan_branch");        
  }
  
  my $job_group = [$fan_branch, $fan_jobs]; 
   
  if(($funnel_branch_code && ! $funnel_jobs) ||
     ($funnel_jobs && ! $funnel_branch_code)){
    $self->throw('Must have both a funnel branch code and funnel job(s) to define a funnel in a job group');    
  }
  elsif($funnel_branch_code){ #implicit that $funnel_jobs are also true
    #check $funnel_branch_code is not ref?
     
    if(! (check_ref($funnel_jobs, 'ARRAY') &&
          scalar(@$funnel_jobs) > 0)){
      $self->throw('Must have at least 1 funnel job when defining a funnel in a job group'); 
    }
   
    my $funnel_branch = $self->_get_branch_number([$funnel_branch_code], $branch_config); 
    $self->helper->debug(1, 'Branching funnel '.scalar(@$funnel_jobs).
      " to branch $funnel_branch(code=$funnel_branch_code)");
   
   
    if($branch_config){#check the funnel is valid wrt fan analyses
      #Funnel branch code, may also just be a branch number and not a logic name
      #So we have to use the funnel_branch and fan branch rather than codes
      my $funnel_name = $branch_config->{$fan_branch}{funnel};
      
      if($branch_config->{$funnel_name}{branch} ne $funnel_branch){
        $self->throw("$funnel_branch_code (expected: ".$branch_config->{$funnel_name}{branch}.") is not a valid funnel analysis for fan branches:\t".
          join(', ', @$fan_branch_codes)."\nPlease check you dataflow configuration");  
      }   
    }
    
    push @$job_group, ($funnel_branch, $funnel_jobs);
  }
  
  $self->{job_groups} ||= [];
  push @{$self->{job_groups}}, $job_group;
  
  return;
}

=head2 dataflow_job_groups

  Arg [1]    : Arrayref - Optional job groups
  Example    : 
  Description: 
  Returntype : 
  Exceptions : 
  Caller     : General
  Status     : At risk

=cut
sub dataflow_job_groups {
  my $self       = shift;
  my $job_groups = $self->{job_groups} || []; #allow no job groups
  
  foreach my $job_group(@$job_groups){
    #we can detect funnel jobs here!
    #Is it possibl to have more than one semaphore from a fan?
    #when will the fan be flown, on the first funnel flow? 
     
     
    while(@$job_group){
      my $branch   = shift @$job_group;
      my $id_array = shift @$job_group;
      
      #let dataflow_output_id do the validation
      
      #This will fail silently if the input_id already exists
      #There is no way to know which of the jobs have failed to flow
      #and it only returns this internal DB ids of the successful jobs
    
      #why are we only getting 1 job in the group here?
      #is this how preprocess flows them?
    
      warn "Dataflowing ".scalar(@$id_array)." jobs to branch $branch";
    
      $self->dataflow_output_id($id_array, $branch);  

    }
  } 

  return;
}

sub dataflow_params {
  my $self         = shift;
  my $optional     = shift;
  my $param_values = {};
  
  if(! $self->param_silent('dataflow_param_names')){
    
    if(! $optional){
      $self->throw('dataflow_param_names are required but not defined.'.
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
  my %param_values = ();#Define as we will deref in the caller
  
  $self->helper->debug(3, "Getting $method_list_name ", $self->param_required($method_list_name)); 
   
  foreach my $param_name(@{$self->param_required($method_list_name)}){
  
    my $param = $self->param_silent($param_name);
    
    if( blessed($param) ){
      
      if(! exists $object_dataflow_methods{ref($param)}){
        $self->throw("Cannot dataflow $method_list_name $param_name as it is an object without a defined dataflow method:\t".
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

sub sam_ref_fai {
  my $self        = shift;
  
  use Bio::EnsEMBL::Funcgen::Hive::RefBuildFileLocator;
  
  my $bwa_index_locator = Bio::EnsEMBL::Funcgen::Hive::RefBuildFileLocator->new;
  
  my $reference_file_root = $self->param('reference_data_root_dir');

  my $samtools_fasta_index_relative = $bwa_index_locator->locate({
    species          => $self->species,
    assembly         => $self->assembly,
    epigenome_gender => $self->param('gender'),
    file_type        => 'samtools_fasta_index',
  });
  
  my $samtools_fasta_index = $reference_file_root . '/' . $samtools_fasta_index_relative;
  
  $self->param('sam_ref_fai', $samtools_fasta_index);
  
  return $samtools_fasta_index;
}

sub get_alignment_path_prefix_by_ResultSet {
  my ($self, $rset, $control) = @_;
  assert_ref($rset, 'Bio::EnsEMBL::Funcgen::ResultSet');

  return $self->alignment_dir($rset, undef, $control).'/'.
    $self->get_file_base_for_ResultSet($rset, $control); 
}

sub get_file_base_for_ResultSet {
  my ($self, $rset, $control) = @_;
  assert_ref($rset, 'Bio::EnsEMBL::Funcgen::ResultSet');
  
  my $ctrl_exp;
  if ($control) {
    $ctrl_exp = $rset->experiment(1);
  }
  return if $control && ! $ctrl_exp;
  
  my @rep_numbers;
  foreach my $isset(@{$rset->get_support}) {
    if( ($control && $isset->is_control)  ||
        ((! $control) && 
         (! $isset->is_control)) ) {
       push @rep_numbers, $isset;
    }
  }
  
  my $replicate_input_subset_string = $self->create_replicate_input_subset_string(@rep_numbers);
  
  my $file_base_for_result_set = get_set_prefix_from_Set($rset, $control)
    . '_'
    . $rset->analysis->logic_name 
    . '_'
    . $replicate_input_subset_string;
  
  return $file_base_for_result_set;
}

sub create_replicate_input_subset_string {

  my $self = shift;
  my @input_subset_replicates = @_;

  my @sorted_replicate_number = sort {
       $a->biological_replicate <=> $b->biological_replicate
    || $a->technical_replicate  <=> $b->technical_replicate
  } @input_subset_replicates;
  
  my @replicate_number_strings;
  foreach my $current_replicate_number (@sorted_replicate_number) {
    push @replicate_number_strings,
        'BR' . $current_replicate_number->biological_replicate
      . 'TR' . $current_replicate_number->technical_replicate;
  }
  my $replicate_number_string = join('_', sort(@replicate_number_strings));
  
  return $replicate_number_string;
}

=head2 signal_input_subsets_have_biological_replicates

  This code is used in two analysis modules and they rely on it returning
  the same result both times.

=cut
sub signal_input_subsets_have_biological_replicates {

  my $self  = shift;
  my $signal_input_subsets = shift;
  
  my $biological_replicates_present = grep {
    $_->biological_replicate > 1;
  } @$signal_input_subsets;

  return $biological_replicates_present;
}

sub create_result_set_name_for_biological_replicate {
  my $self  = shift;
  my $param = shift;
  
  my $parent_result_set_name = $param->{parent_result_set_name};
  my $replicate_number       = $param->{replicate_number};
  
  return $self->_create_result_set_name_for_replicate(
    $parent_result_set_name, 
    'BR', 
    $replicate_number
  );
}

sub create_result_set_name_for_technical_replicate {
  my $self  = shift;
  my $param = shift;
  
  my $parent_result_set_name = $param->{parent_result_set_name};
  my $replicate_number       = $param->{replicate_number};
  
  return $self->_create_result_set_name_for_replicate(
    $parent_result_set_name, 
    'TR', 
    $replicate_number
  );
}

sub _create_result_set_name_for_replicate {
  my $self  = shift;
  
  my $parent_result_set_name = shift;
  my $replicate_type         = shift;
  my $replicate_number       = shift;
  
  return $parent_result_set_name . '_' . $replicate_type . $replicate_number;
}

sub get_alignment_files_by_ResultSet_formats {
  my ($self, $result_set, $control) = @_;
   my ($path, $align_files);
 
  if (!$control) {
    $path = $result_set->dbfile_path;
  }
  
  if ($control) {
  
    my $control_experiment = $result_set->experiment(1);
    
#     if ($experiment->is_control) {
#       $self->throw("This experiment (". $experiment->name .") from the result set " . $result_set->name . " is already the control. It does not have a control.");
#     }
#     my $control_experiment = $experiment->get_control;

    my $control_result_set_array_ref = $result_set->adaptor->fetch_all_by_Experiment($control_experiment);
    
    die("Couldn't find control result set!") unless($control_result_set_array_ref);
    die("Couldn't find unique control result set!") unless(@$control_result_set_array_ref>1);
    
    my $control_result_set = $control_result_set_array_ref->[0];
    
#     use Data::Dumper;
#     die( Dumper($control_result_set->[0]) );
#     die;
    
    $path = $control_result_set->dbfile_path;
  }
  
  return $path;

}

=head2 is_idr_FeatureType

  Everything is an idr feature type unless it is a broad peak feature type.
  
  If no_idr has been set, nothing is an idr feature type.

=cut
sub is_idr_FeatureType {
  my $self  = shift;
  my $ftype = shift;
  assert_ref($ftype, 'Bio::EnsEMBL::Funcgen::FeatureType');
  $ftype = $ftype->name;
  
  # This ! grep is fine, although it returns an empty string instead of 0
  # as oppose to 1, when nothing is returned from grep
  #
  return $self->param_silent('no_idr') ? 0 : 
          ! grep(/^${ftype}$/, @{$self->param_required('broad_peak_feature_types')});
}

sub is_idr_ResultSet {
  my $self = shift,
  my $rset = shift;
  assert_ref($rset, 'Bio::EnsEMBL::Funcgen::ResultSet');
  my $is_idr_rset = 0;  
  my @sig_reps    = grep { ! $_->is_control } @{$rset->get_support}; 
    
  if($self->is_idr_FeatureType($rset->feature_type) &&
     (scalar(@sig_reps) > 1) ){
    $is_idr_rset = 1;     
  }

  return $is_idr_rset;  
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
  $self->throw(@_);

} ## end throw_no_retry  

1;
