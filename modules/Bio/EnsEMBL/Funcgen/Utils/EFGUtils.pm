=head1 LICENSE

  Copyright (c) 1999-2013 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=head1 NAME

Bio::EnsEMBL::Funcgen::Utils::EFGUtils

=head1 DESCRIPTION

This module collates a variety of miscellaneous methods.


=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(list your required methods here);



=cut


# No API/Object based methods in here

###############################################################################

package Bio::EnsEMBL::Funcgen::Utils::EFGUtils;

require Exporter;

@ISA = qw(Exporter);

@EXPORT_OK = qw(
  add_external_db
  add_hive_url_to_meta      
  are_valid
  assert_ref_do     
  assert_refs               
  assert_refs_do  
  backup_file
  create_Storable_clone 
  generate_slices_from_names
  get_current_regulatory_input_names          
  get_date
  get_feature_file
  get_file_format           
  get_month_number
  get_my_scalar_name        
  get_feature_file  
  is_bed            
  is_gzipped  
  is_sam        
  mean              
  median        
  open_file   
  parse_DB_url      
  run_system_cmd 
  scalars_to_objects
  set_attributes_by_my_scalar_names
  species_chr_num
  species_name
  strip_param_args  
  strip_param_flags
  url_from_DB_params
  validate_path
  );

#These after @ISA and @EXPORT_OK to avoid warnings               
use warnings;
use strict;
               
use Bio::EnsEMBL::Utils::Exception qw( throw      );
use Bio::EnsEMBL::Utils::Scalar    qw( assert_ref );
use Scalar::Util                   qw( blessed    );
use File::Path                     qw( make_path  );
use File::Basename                 qw( dirname fileparse );
use File::Spec;
use Time::Local;
use FileHandle;
use Carp;


#Handy function to add an external_db entry when needed
sub add_external_db{
  my ($efg_db, $db_name,$db_release,$db_display_name) = @_;
  my $sql = "select external_db_id from external_db where db_name='$db_name' and db_release='$db_release'";
  my ($db_id) =  $efg_db->dbc->db_handle->selectrow_array($sql);
  if($db_id){
    warn "External DB $db_name $db_release already exists in db with db_id $db_id\n";
  } else {
    #TODO check if it there was a failure
    $efg_db->dbc->do("insert into external_db (db_name, db_release, status, dbprimary_acc_linkable, priority, db_display_name, type) values('$db_name', '$db_release', 'KNOWNXREF', '1', '5', '$db_display_name', 'MISC')");
  }
}


#Move to PipelineUtils?
 
sub add_hive_url_to_meta {
  my ($url, $db) = @_;  
 
  #VALIDATE ARGS
  parse_DB_url($url);
  assert_ref($db, 'Bio::EnsEMBL::DBSQL::DBAdaptor', 'db');#restrict this to Funcgen?
  
  my $mc         = $db->get_MetaContainer;
  my $meta_key = 'hive_url';
  
  my $meta_value = $mc->single_value_by_key($meta_key);

  if( ! defined $meta_value){ #Store the new meta entry    
    #Add key via API to store with appropriate species_id
    eval { $mc->store_key_value($meta_key, $url); };
    
    if($@){
      throw("Failed to store hive url meta entry:\t$meta_key\t$url\n$@");  
    }
  }
  elsif($meta_value ne $url){
    throw("Could not set hive_url meta entry:\t".$url."\nAs ".$db->dbc->dbname.
      " is currently locked to a different hive DB:\t$meta_value\n".
      'Please use that hive DB or drop that pipeline, or remove the meta entry');
  }
  #else meta values match and all is good.
  
  return;
}


=head2 are_valid

  Arg [1]    : String - Namespace of object class
  Arg [2]    : ARRAYREF of objects to be tested
  Arg [3]    : String (optional) - Method name to call on each object.
  Arg [4]    : Bio::EnsEMBL::Funcen::DBSQL::DBAdaptor (optional)
  Example    : $db->are_valid('Bio::EnsEMBL::Funcgen::ResultSet', \@rsets, undef, $db);
  DESCRIPTION: Validate objects namespaces and optionally (if a DBADaptor is passed)
               whether they are stored. If method name arg is defined, this will be called
               on each object to populate the return value array.
  Returntype : ARRAYREF - contents defined by optional method name arg
  Exceptions : Throws if object list is not an ARRAY with at least one element
               Throws if an object does not inherit from the given class
  Caller     : General
  Status     : At risk

=cut

#add method params here?

sub are_valid {
 my ($class, $obj_list, $method_name, $db) = @_;

  assert_ref($obj_list, 'ARRAY', 'object list');
  
  if(scalar(@$obj_list) == 0){
   throw('Objects Arrayref is empty'); 
  }
  
  my @return_vals;
 
  foreach my $obj (@$obj_list) {
    
    if($db){
      $db->is_stored_and_valid($class, $obj);
    }
    
    if(! $method_name){
      assert_ref($obj, $class, 'object');
    }
    else{
      push @return_vals, assert_ref_do($obj, $class, $method_name, 'object');
    }
  }

  return \@return_vals;
}
 


=head2 assert_ref_do

  Arg [1]     : The reference to check
  Arg [2]     : The method we want to call
  Arg [3]     : The attribute name you are asserting; not required but allows
                for more useful error messages to be generated. Defaults to
                C<-Unknown->.
  Description : A subroutine which checks to see if the given object/ref is
                implements the given method. Will throw exceptions.
                
  Returntype  : Boolean; true if we managed to get to the return
  Example     : assert_ref_can($gene, 'dbID');
  Exceptions  : If the reference is not defined, if the object does not
                implement the given method and if no method was given to check
  Status      : Stable

=cut

#Ideally want to make assert_ref support can and do
#then we can make all these wrapper to that method passing a $do_can arg

sub assert_ref_do {
  my ($ref, $expected, $method, $attribute_name) = @_;

  #Can't do this as we expect a method return value
  #return 1 unless $ASSERTIONS;
  my $value;
  
  $attribute_name ||= '-Unknown-';
  throw('No method given') if ! defined $method; 
  throw('No expected type given') if ! defined $expected;
  throw "The given reference $attribute_name is not defined" unless defined $ref;
  throw "The given reference $attribute_name is not blessed" unless blessed($ref); 
  throw("${attribute_name}'s type '${ref}' is not an ISA of '${expected}'") if ! $ref->isa($expected);

  if(! $ref->can($method)) {
    my $str_ref = ref($ref);
    throw sprintf(q{The given blessed reference '%s' for attribute '%s' does not implement the method '%s'}, $str_ref, $attribute_name, $method);
  }
  
  
  return $ref->$method;
}


#Return type isn't actually boolean as it is never 0

sub assert_refs {
  my ($refs, $expected, $attribute_name) = @_;
  
  my $array_attr_name = (defined $attribute_name) ? 
    $attribute_name.' array' : undef;
    
  assert_ref($refs, 'ARRAY', $array_attr_name); 
  map { assert_ref($_, $expected, $attribute_name) } @$refs;
   
  return 1; 
}


=head2 assert_refs_do

  Arg [1]     : An Arrayref of references to check
  Arg [2]     : The type we expect
  Arg [3]     : The attribute name you are asserting; not required but allows
                for more useful error messages to be generated. Defaults to
                C<-Unknown->.
  Description : A subroutine which checks to see if the given objects/refs are
                what you expect. This behaves in an identical manner as
                C<assert_ref> does works on an array ref of references
                
                You can turn assertions off by using the global variable
                $Bio::EnsEMBL::Utils::Scalar::ASSERTIONS = 0
  Returntype  : Boolean; true if we managed to get to the return
  Example     : assert_refs([[],[],[]], 'ARRAY');
  Exceptions  : Throws is references argument is not an Arrayref, also
                if the expected type was not set and if the given reference
                was not assignable to the expected value.
  Status      : Stable

=cut


sub assert_refs_do {
  my ($refs, $expected, $method, $attribute_name) = @_;
  
  my $array_attr_name = (defined $attribute_name) ? 
    $attribute_name.' array' : undef;
    
  assert_ref($refs, 'ARRAY', $array_attr_name); 
  my @values;
  map { push @values, assert_ref_do($_, $expected, $method, $attribute_name) } @$refs;
   
  return \@values; 
}


#can't do check_ref_do as this would make the return type ambiguous


sub backup_file{
  my $file_path = shift;

  throw("Must define a file path to backup") if(! $file_path);

  if (-f $file_path) {
    #$self->log("Backing up:\t$file_path");
    system ("mv ${file_path} ${file_path}.".`date '+%T'`) == 0 || return 0;
  }

  return 1;

}



=head2

  Name       : create_Storable_clone
  Arg [1]    : Bio::EnsEMBL::Funcgen::Storable to clone
  Arg [2]    : Hash containing list of object parameters linked to
               the clone and to be reset
  Example    :
  Description: Blesses an object and replaces all linked objects with the ones
               passed with the $params_hash
  Returntype : cloned object
  Exceptions : Throws if object is not a Storable
  Caller     : general
  Status     : At risk - not tested

=cut

#Currently not exposing the no db reset flag from reset_relational_attributes
#as we would never want to do this here

sub create_Storable_clone {
  my ($obj, $params_hash) = @_;
  
  assert_ref($obj, 'Bio::EnsEMBL::Funcgen::Storable', 'object');
  
  if(! $obj->can('reset_relational_attributes')){
     throw('You must pass a Bio::EnsEMBL::Funcgen::Storable which can '.
      'call the reset_relational_attributes method');
  }

  my $clone = bless({%{$obj}}, ref($obj));
  $clone->reset_relational_attributes($params_hash);
  return $clone;
}


#Do this as File::Basename doesn't quite get there.

sub file_suffix_parse{
  my $filepath = $_[0];
  
  my $file_name                  = fileparse($filepath);
  (my $file_prefix = $file_name) =~ s/\.[a-zA-Z0-9]$//;
  (my $suffix      = $file_name) =~ s/$file_prefix\.//;
    
  return ($file_prefix, $suffix);
}


#Generates slices from names or optionally alll default top level nonref
#slice ref args can be array ref (inc empty) or undef

sub generate_slices_from_names{
  my ($slice_adaptor, $slice_names, $skip_slices, $level, $non_ref, $inc_dups, $assembly) = @_;
  my (@slices, $slice, $sr_name, $have_slice_names, $have_skip_slices);

  #Test if $assembly is old?

  #Validate array ref skip/slice name args

  if(defined $slice_names){
    if(ref($slice_names) ne 'ARRAY'){
      throw('Slice names argument must be an ARRAYREF');
    }

    $have_slice_names = 1 if @$slice_names;
  }

  if(defined $skip_slices){
    if(ref($skip_slices) ne 'ARRAY'){
      throw('Skip slices argument must be an ARRAYREF');
    }

    $have_skip_slices = 1 if @$skip_slices;
  }

  #Generate slices

  if($have_slice_names){

    foreach my $name(@$slice_names){
      $slice = $slice_adaptor->fetch_by_region(undef, $name, undef, undef, undef, $assembly);

      #WHy is this failing for hap regions?

      if(! $slice){

        #Need to eval this as it will break with incorrect formating

        eval { $slice = $slice_adaptor->fetch_by_name($name) };

        if(! $slice){
          throw("Could not fetch slice by region or name:\t".$name);
        }
      }

      $sr_name = $slice->seq_region_name;

      next if(grep/^${sr_name}$/, @$skip_slices);
      push @slices, $slice;
    }
  }
  elsif($level){

    #Can't guarantee what level will have an assembly version
    #Can only throw if we have set toplevel and assembly
    #as this will alway default to the current assembly
    #if($assembly &&
    #   ($level eq 'toplevel') ){
    #  throw('You cannot specify an assembly version with the toplevel coordinate system');
    #}
    #Let core API handle this?


    my @tmp_slices = @{$slice_adaptor->fetch_all($level, $assembly, $non_ref, $inc_dups)};

    if($have_skip_slices){

      foreach $slice(@tmp_slices){
        $sr_name = $slice->seq_region_name;
        push @slices, $slice if ! grep/^${sr_name}$/, @$skip_slices;
      }
    }
    else{
      @slices = @tmp_slices;
    }
  }
  else{
    throw('You must either pass an arrayref of slice names or specify the toplevel flag');
  }


  if(! @slices){
    throw("You have specified slice_names and skip_slices paramters which have generated no slices.\nslice_names:\t".join(' ',@$slice_names)."\nskip_slices:\t".join(' ', @$skip_slices));
  }

  return \@slices;
}


# Tracking DB methods
# Move to DBAdaptor? Can we add this as a separate package in the same module?

sub get_current_regulatory_input_names{
  my ($tdb, $efg_db, $focus) = @_;

  warn "Move get_current_regulatory_input_names to the TrackingAdaptor";

  #Validate is production?
  my $sql;



  if($focus){
    $focus = 'Focus';
    $sql   = 'SELECT efgdb_set_name from dataset where is_focus=true and is_current=true and species="'.$efg_db->species.'"';
  }
  else{
    $focus = 'Non-focus';
    #0 rather than false so we don't get NULLs
    $sql = 'SELECT efgdb_set_name from dataset where is_focus=0 and is_current=true and species="'.$efg_db->species.'"';
  }



  #Currently efgdb_set_name can either be data_set or feature_set name!
  #Need to standardise this

  my @prd_names = @{$tdb->db_handle->selectcol_arrayref($sql)};
  my @names;
  my @failed_sets;

  foreach my $prd_name(@prd_names){

    $sql = "SELECT name from feature_set where name like '${prd_name}%'";
    my @tmp_names =  @{$efg_db->dbc->db_handle->selectcol_arrayref($sql)};

    #This is causing problems with multiple feature sets with differing analyses

    #Do this via InputSets(using query extension?) instead of using like?

    #This is very hacky right now to get it to work
    #Need to standardise and review tracking db data.

    if(scalar(@tmp_names) > 1){

      $sql = "SELECT name from feature_set where name ='${prd_name}_ccat_histone'";
      @tmp_names =  @{$efg_db->dbc->db_handle->selectcol_arrayref($sql)};

      if(scalar(@tmp_names) == 1){
        push @names, $tmp_names[0];
      }else{
        push @failed_sets, $prd_name;
      }
    }
    elsif(scalar(@tmp_names) == 0){
      push @failed_sets, $prd_name;
    }
    else{
      push @names, $tmp_names[0];
    }
  }

  if(@failed_sets){
    throw("Failed to find unique $focus FeatureSets for production dataset names:\n\t".
          join("\n\t", @failed_sets)."\n");
  }

  return @names;
}


sub get_date{
	my ($format, $file) = @_;

	my ($time, $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst);


	throw("File does not exist or is not a regular file:\t$file") if $file && ! -f $file;


	($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = (defined $file) ?
	  localtime((stat($file))[9]) : localtime();

	#print "	($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst)\n";

	if((! defined $format && ! defined $file) || $format eq "date"){
		$time = ($year+1900)."-".$mday."-".($mon+1);
	}
	elsif($format eq "time"){#not working!
		$time = "${hour}:${min}:${sec}";
	}
	elsif($format eq "timedate"){#
	  $time = localtime();
	}
	else{#add mysql formats here, datetime etc...
		croak("get_date does not handle format:\t$format");
	}

	return $time;
}



#This handles g/unzipping and conversion from bam > sam > bed
#This also assumes that we only ever want to convert in this direction
#i.e. assumes bam /sam will always exist if we have bed.

#sam params contains:
#ref_fai         => file_path
#filter_from_bam => 1
#could also support:
#include_MT   => 1,
#include_dups => 1,
#ignore header mismatch? (could do thi swith levels, which ignore supersets in fai?

#Currently hardoded for samse files name for sam and bed
#todo _validate_sam_params

sub get_feature_file {
  my ($path, $sam_params) = @_;  
  #warn "$path $bam_params";  
  my $required_suffix = $path;
  $required_suffix =~ s/.*\.//;
  my $last_suffix     = $required_suffix;
  my $tmp_path        = $path;
  my @suffixes;
 
  if(! defined $required_suffix){
    throw("Unable to identify file format/suffix from:\t$path"); 
  }
  
  if($sam_params && (ref($sam_params) ne 'HASH') ){
    throw('sam params argument must be a Hashref');
  }
  #else validate keys here
  
  if($sam_params->{filter_from_bam}){
    
    #Do we expect a full/unfiltered file to be present?
    #or do we just use the existing file
    
    if($required_suffix eq 'bam'){
      throw('get_feature_file does not yet support forcing filtering from bam if the target format is also bam');
    }
  
    @suffixes = ('bam');
  }
  else{ 
    @suffixes = ($required_suffix);
  
    if($required_suffix eq 'sam'){
      push @suffixes, 'bam';
    }
    elsif($required_suffix eq 'bed'){
      push @suffixes,  ('sam', 'bam');
    }
    elsif($required_suffix ne 'bam'){
      throw("$required_suffix is not a recognised feature file format. Must be bam, sam or bed");  
    }
  }

  #Try all variants of suffixes and gzipped files
  my ($found_suffix, $found_path);
 
  foreach my $suffix(@suffixes){
    #my $sub_suffix = ($suffix eq 'bam') ? "\.samse\.$last_suffix" : '\.'.$last_suffix;
    #Shouldn't our bam files also have samse in them? YES! This may have been chaged recently, but should be re-implemented
    #or should it? This will add the required for a samse/sampe flag somewhere?
    #omission will mean we won't know from the data file named, whether it was a single end or a paired end alignment
    #actually a filtered file based on quality score can contain both se and pe reads
    $tmp_path =~ s/\.$last_suffix$/\.$suffix/; # Define path to be tried   
    $found_suffix = $suffix;
    
    if(! -f $tmp_path){                      # Try gzipped file
         
      if(-f $tmp_path.".gz"){                # Gotcha!  
        $found_path = $tmp_path;
      
        if($suffix eq $required_suffix){
          run_system_cmd("gunzip ${tmp_path}.gz");
        }
        else{
          $found_path .= ".gz";
        }         
        last;
      }
    
      $last_suffix = $suffix;      
      #This is the only point where we actually iterate       
    }
    else{                                    # Gotcha!
      $found_path = $tmp_path;
      last; 
    }                           
  
    #Should never get to this point
  }
  
  
  
  if(! defined $found_path){
    throw("Could not find feature file:\t$path\nOr useable ".join(' OR ', @suffixes).' (gz) files');   
  }
  
  #Need to handle rezipping of files!
  
  #("$required_suffix ne $found_suffix found file $found_path");
  
  if($required_suffix ne $found_suffix){
    
    #Need to do a FindBin or similar here to get $EFG_SRC rather
    #than relying on env var
    
    
    if($found_suffix eq 'bam'){ #Do bam2sam
      $found_path   = sort_and_filter_sam($found_path, $sam_params->{ref_fai});
      $found_suffix = 'sam';
    }
    
     #Use bam to bed from BedTools instead?
     #This is another pre-req!
    
    if($required_suffix ne $found_suffix){ #required suffix has to be bed and suffix will be sam
    
      #if($suffix eq 'sam'){         # sam2bed
        
      #sam2bed handles gzip files
      run_system_cmd($ENV{EFG_SRC}."/scripts/miscellaneous/sam2bed.pl -uncompressed -1_based -files $found_path");
      
      if($found_suffix eq 'bam'){    #Delete intermediate sam file
        run_system_cmd("rm -f $found_path");
      }
       
      $found_path =~ s/\.sam(\.gz)*$/.bed/;    
    }   
  }
  
  return $found_path;
}


sub get_file_format{
  my $file = shift;

  my $format = &is_bed($file);

  if(! $format){
    $format =  &is_sam($file);

    #Add more tests here
  }


  return $format;
}



sub get_month_number{
  my($mon) = @_;
  my %month_nos =(
          "jan", "01",
          "feb", "02",
          "mar", "03",
          "apr", "04",
          "may", "05",
          "jun", "06",
          "jul", "07",
          "aug", "08",
          "sep", "09",
          "oct", "10",
          "nov", "11",
          "dec", "12",
         );
  return (exists $month_nos{lc($mon)}) ? $month_nos{lc($mon)} : undef;
}


sub gunzip_file {
  my ($self, $filepath) = @_;
  my $was_gzipped = 0;
  
  if( is_gzipped($filepath) ){
    system("gunzip $filepath") && throw("Failed to gunzip file:\n$filepath\n$@"); 
    $filepath =~ s/\.gz$//; 
    $was_gzipped = 0;
  }

  return($filepath, $was_gzipped);
}


sub is_bed {
  my $file = shift;

  #Use open_file here!

  if(&is_gzipped($file, 1)){

    open(FILE, "zcat $file 2>&1 |") or throw("Can't open file via zcat:\t$file");
  }
  else{
    open(FILE, $file) or throw("Can't open file:\t$file");
  }

  my @line;
  #$verbose =1;


  while (<FILE>) {
    chomp;
    @line = split("\t", $_);
    last;
  }
  close FILE;

  if (scalar @line < 6) {
    #warn "Infile '$file' does not have 6 or more columns. We expect bed format:\t".
    #  "CHROM START END NAME SCORE STRAND.\n";
    return 0;
    #} elsif ($line[0] !~ m/^((chr)?[MTXYNT_\d]+)$/) {
    #    warn ("1st column must contain name of seq_region (e.g. chr1 or 1) in '$file'");
    #    return 0;
    #Commented this out for now due to HSCHR_RANDOM seqs
    #How does the webcode handle this?
  }
  elsif ($line[1] !~ m/^\d+$/ && $line[2] =~ m/^\d+$/) {
    warn "2nd and 3rd column must contain start and end respectively in '$file'\n";
    return 0;
  }
  elsif ($line[5] !~ m/^[+-\.]$/) {
    warn "6th column must define strand (either +, - or .) in '$file'\n";
    return 0;
  }

  return 'bed';
}


sub is_gzipped {
  my ($file, $fail_if_compressed) = @_;

  throw ("File does not exist:\t$file") if ! -e $file;

  open(FILE, "file -L $file |")
    or throw("Can't execute command 'file' on '$file'");
  my $file_info = <FILE>;
  close FILE;

  my $gzip = 0;

  if($file_info =~ /compressed data/){

    if($file_info =~ /gzip/){
      $gzip = 1;
    }
    else{
      throw("File is compressed but not with gzip, please unzip or gzip:\t$file_info");
    }
  }

  return $gzip;
}



#Change these to also return the gz status

sub is_sam{
  my $file = shift;

  #warn "Only checking file suffix for is_sam";
  #Could check for header here altho this is not mandatory!
  #Can we use web format guessing code?

  my $gz = (&is_gzipped($file, 1)) ? '.gz' : '';

  return ($file =~ /.sam${gz}/) ? 'sam' : 0;
}

#need is bam here too!


sub mean{
  my $scores = shift;

  my $total = 0;

  map $total+= $_, @$scores;
  my $mean = $total/(scalar(@$scores));

  return $mean;

}


#Sort should always be done in the caller if required

sub median{
  my ($scores, $sort) = shift;

  return undef if (! @$scores);


  my ($median);
  my $count = scalar(@$scores);
  my $index = $count-1;
  #need to deal with lines with no results!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #deal with one score fastest
  return  $scores->[0] if ($count == 1);


  if($sort){
    #This is going to sort the reference here, so will affect
    #The array in the caller
    #We need to deref to avoid this
  }

  #taken from Statistics::Descriptive
  #remeber we're dealing with size starting with 1 but indices starting at 0

  if ($count % 2) { #odd number of scores
    $median = $scores->[($index+1)/2];
  }
  else { #even, get mean of flanks
    $median = ($scores->[($index)/2] + $scores->[($index/2)+1] ) / 2;
  }


  return $median;
}


sub open_file{
  my ($file, $operator, $file_permissions) = @_;

  $operator ||= '<';

  if ($operator !~ /%/) {
    $operator = "$operator $file";
  } else {
    #We have some piping to do
    $operator = sprintf($operator, $file);
  }

  #Get dir here and create if not exists
  my $dir = dirname($file);

  my $mkpath_opts = {verbose => 1};
  $mkpath_opts->{mode} = $file_permissions if defined $file_permissions;


  if(! -d $dir){
    eval { make_path($dir, $mkpath_opts) };
    
    if($@){
      throw("Failed to make_path:\t$dir\n$@");
    }
  }
 
  #my $fh;
  
  #Test exists and is readable/writable first to avoid zombie?
  
  my $fh = new FileHandle "$operator";
  #warn "operator is $operator";

  #open(my $fh, $operator) or die("Failed to open:\t$operator"); 

  

  #This does not catch errors when piping then redirecting
  #as the progam will have opened sucessfully
  #i.e. gzip -c | to file we don't have permission to create!
  #This simply carries on with the reset of the script
  #won't be caught until we close the file descriptor
  #Which is a little late, and we don't catch error codes on close yet
 
  if (! defined $fh) {
    croak("Failed to open $operator");
  }

  #it maybe better to use PerlIO libs here e.g. PerlIO::gzip


  #Have to chmod here as umask will over-ride permissions passed to FileHandle
  if (defined $file_permissions) {

    #Catch non-numeric here as chmod still returns true
    if ($file_permissions =~ /[^0-9]/) {
      croak("Failed to change $file permissions using:\t$file_permissions");
    }

    #chmod requires a literal octal number e.g. 0775 not '0775'
    #should catch numbers as strings here, but perl makes this very hard to test
    #Can't even system this as if we build the cmd line with an octal it will be converted to a decimal
    #There is still no way of testing for non-octal number or string
    #eval/sprintf will also not fail if there are non-octal digits i.e. 1999

    #eval will treat octal number and string as true octal number
    #else will pass non-octal string/number which we can't catch
    chmod(eval($file_permissions), $file);
  }

  return $fh;
}



sub parse_DB_url {
  my $url = $_[0];
  
  my ($dbtype, $user, $pass, $host, $dbname, $port);

  if($url =~ /(.*):\/\/(.*):(.*)@(.*):(.*)\/(.*)/){
    ($dbtype, $user, $pass, $host, $port, $dbname) = ($1, $2, $3, $4, $5, $6);
     
    if(! ($user && $host && $dbname)){
      die("The DB url must contain at least user, host and dbname elements i.e.\n\t".
        'DBTYPE://USER:[PASS]@DBHOST:[PORT]/DBNAME');
    }
  }
  else{
    throw("Unrecognized DB url format:\t$url\nShould be:\n\t".
      'DBTYPE://USER:[PASS]@DBHOST:[PORT]/DBNAME'); 
  } 

  return [$user, $pass, $host, $dbname, $port, $dbtype];
}



################################################################################

=head2 run_system_cmd

 Description : Method to control the execution of the standard system() command

 ReturnType  : none

 Example     : $Helper->debug(2,"dir=$dir file=$file");

 Exceptions  : throws exception if system command returns none zero

=cut

################################################################################

#Allow no_exit as some programs(tab2mage) give successful non-zero exit codes!

sub run_system_cmd{
  my ($command, $no_exit) = @_;
  

  my $redirect = '';

  #$self->debug(3, "system($command)");

  # decide where the command line output should be redirected

  #This should account for redirects

  #if ($self->{_debug_level} >= 3){

  #  if (defined $self->{_debug_file}){
  #    $redirect = " >>".$self->{_debug_file}." 2>&1";
  #  }
  #  else{
  #    $redirect = "";
  #  }
  #}
  #else{
    #$redirect = " > /dev/null 2>&1";
  #}

  # execute the passed system command
  my $status = system("$command $redirect");
  my $exit_code = $status >> 8;


  if ($status == -1) {
    warn "Failed to execute: $!\n";
  }
  elsif ($status & 127) {
    warn sprintf("Child died with signal %d, %s coredump\nError:\t$!",($status & 127),($status & 128) ? 'with' : 'without');
  }
  elsif($status != 0) {
    warn sprintf("Child exited with value %d\nError:\t$!\n", $exit_code); #get the true exit code
  }

  #We're not catching error message here?

  if ($exit_code != 0){

    if (! $no_exit){
      throw("System command failed:\t$command");
    }
    else{
      warn("System command returned non-zero exit code:\t$command");
    }
  }

  return $exit_code;
}


=head2 scalars_to_objects

  Arg [1]    : Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor
  Arg [2]    : String - class name of object to retrieve
  Arg [3]    : String - method name to use
  Arg [4]    : Arrayref - Scalar arguments to use iteratively with the fetch method
  Example    : my @cell_types = @{scalars_to_object($db, 'CellType',
                                                    'fetch_by_name',
                                                    [ qw ( GM06990 HUVEC H1ESC ) ])};
  Description: Convenience method to fetch objects iteratively based on an
               array of scalar arguments. Useful for processing command line arguments.
  Returntype : Arrayref of Objects
  Exceptions : Throws if arguments not specified or valid
               Throws if cannot get relevant adaptor based on the class name
               Throws if Object Adaptor cannot call fetch method
               Throws if Object not retrieved
  Caller     : General
  Status     : At risk

=cut

#We could be clever here and support incorporation of this sub
#into a class by assuming that $db might also be $self if called as a method.
# i.e. could test if can('db') and use that?

sub scalars_to_objects{
  my ($db, $class_name, $fetch_method, $scalars) = @_;
  assert_ref($db, 'Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor', 'db');
  
  if(! ((defined $scalars && (ref($scalars) eq 'ARRAY')) &&
         defined $class_name &&
         defined $fetch_method) ){
    throw('You need to specify a class name and a fetch method and an Arrayref of scalar method arguments');
  }

  my @objs;
  my $adaptor_method = 'get_'.$class_name.'Adaptor';

  my $adaptor = $db->$adaptor_method;

  if(! ((defined $adaptor) &&
        $adaptor->can($fetch_method))){
    throw("Could not $adaptor_method or ${class_name}Adaptor cannot call method $fetch_method");
  }

  foreach my $str(@$scalars){

    my $obj = $adaptor->$fetch_method($str);

    if(! defined $obj){
      throw("Could not fetch object using ${class_name}Adaptor->${fetch_method}('$str')");
    }

    push @objs, $obj;
  }

  return \@objs;
}



#enable sam2bed here too?


#sorts, filters duplications and MT mappings
#todo make filtering optional/configurable
#header now optional as default will be to include it with everything for safety
#How will this behave with an absent in file header and not fai supplied?

#would be nice to add a DESTROY method to remove any tmp sorted files which may persist
#after an ungraceful exit
#THese can be added to a global $main::files_to_delete array
#which should then also be undef'd in DESTROY so they don't persisnt to another instance

#how are we maintaining the original unfiltered file? See pipeline dev notes about this
#and the 'bam_filtered' config

sub sort_and_filter_sam {
  my ($sam_bam_path, $fasta_fai, $out_format, $out_file) = @_;
      
  validate_path($sam_bam_path);
  my $fasta_fai_opt = '';
  
  if (defined $fasta_fai) {
    validate_path($fasta_fai);
    $fasta_fai_opt = " -t $fasta_fai ";
    #-t spec is now optional, as not integrating the header into each file
    #is risky as it could get regenerated and be mismatched
    #which may produce unpredicatable behaviour or faults if the header
    #doesn't match the data. 
  }
 
  my $out_flag   = ''; #sam default
  $out_format  ||= 'sam';
   
  if($out_format eq 'bam'){
    $out_flag = 'b'; 
  }
  elsif($out_format ne 'sam'){
    throw("$out_format is not a valid samtools output format"); 
  }
   
   my $in_format = 'sam';
   my $in_flag   = 'S';
  
  if($sam_bam_path =~ /\.bam$/){     # bam
    #bam should always have header? Is this true?
    #$cmd = "samtools view -h $sam_bam_path | ";
    ($out_file ||= $sam_bam_path) =~ s/\.bam$/.filtered.${out_format}/;
    $in_format = 'bam';
    $in_flag   = '';
  } 
  elsif($sam_bam_path =~ /\.sam\.gz$/){ # zipped sam
    #$cmd = "gzip -dc $sam_bam_path | ";
    $sam_bam_path =~ s/\.gz$//;
    ($out_file ||= $sam_bam_path) =~ s/\.sam\.gz$/.filtered.${out_format}/;   
  }
  elsif($sam_bam_path =~ /\.sam$/){  # unzipped sam
    ($out_file ||= $sam_bam_path) =~ s/\.sam$/.filtered.${out_format}/;   
    #$cmd = "cat $sam_bam_path | ";
  }
  else{
    throw("Unrecognised sam/bam file:\t".$sam_bam_path); 
  }

  #Could do all of this with samtools view?
  #Will fail if header is absent and $fasta_fai not specified 
  #$fasta_fai is ignored if infile header present
  #Doing it like this would integrate the fai into the output file, which is probably what we want?
  #This would also catch the absent header in the first command rather than further down the pipe
  #chain which will not becaught gracefully

  #This can result in mismatched headers, as it does seem like the fai file is used rather than the in file header
  #here, at least for bams
  
  #We need to test whether it has a header or not before we do this!
  my ($fai_header, $in_file_header);
  
  eval { $in_file_header =`samtools view -H $sam_bam_path` };
  
  #$! not $@ here which will be null string
  
  if($!){
    
    if(! defined $fasta_fai){
      throw("Could not find an in file header or a sam_ref_fai for:\t$sam_bam_path$!");
    }
  }
  elsif(defined $fasta_fai){
    #Now what do we do?
    #will -H now return the faig header, or the infile header?
    #I think the fai header (at least for bam) 
    
    $fai_header = `samtools view -H $fasta_fai_opt $sam_bam_path`; 
    
    if($!){
      throw("Failed to identify view valid header from:\t$fasta_fai\n$!"); 
    }
    
    #Just make sure they are the same
    if($fai_header ne $in_file_header){
      warn("Found mismatched infile and fai headers for:\n\t$sam_bam_path\n\t$fasta_fai");
      #should this throw or just default to in file?  
      #shoudl probably give the options to reheader and sub out this whole compare_sam_header thing
      #reheader mode, would that all old @SQ are present in new fai header?
      #can we get the header form the fai by simply providing an empty file?
      #cat '' | samtools view -HSt fasta.fai -
      
    }
    
    #Just undef this for simplicity
    $fasta_fai_opt = '';
  }
  


  my $cmd = "samtools view -h${in_flag} $fasta_fai_opt $sam_bam_path -F 4 | ". # Incorporate header into file and filter unmapped reads
    " grep -vE '^[^[:space:]]+[[:blank:]][^[:space:]]+[[:blank:]][^[:space:]]+\:[^[:space:]]+\:MT\:' ". #Filter MTs
    " | grep -v '^MT' | grep -v '^chrM' | ";                                                            #Filter more MTs


  #Remove unmapped reads...  
  #-u uncompressed bam (as we are piping)
  #-S SAM input
  #-t  header file (could omit this if it is integrated into the sam/bam?)
  #-F Skip alignments with bit set in flag (4 = unaligned)
  #- (dash) here is placeholder for STDIN (as samtools doesn't fully support piping).
  #This is interpreted by bash but probably better to specify /dev/stdin?
  #for clarity and as some programs can treat is as STDOUT or anything else?
  #-b output is bam  
  (my $sorted_bam_prefix = $sam_bam_path) =~ s/\.$in_format//;
  $sorted_bam_prefix .= ".sorted_$$";
  
  #$cmd .= "samtools view -uShb $fasta_fai -F 4 - | ". #Now done above
  $cmd .= "samtools view -uShb - | ".               #simply convert to bam using infile header
    "samtools sort - $sorted_bam_prefix";      #and finally sort
  #-m 2000000 (don't use 2G here,a s G is just ignored and 2k is used, resulting in millions of tmp files)
  #do we need an -m spec here for speed? Probably better to throttle this so jobs are predictable on the farm
  
  #We could also test the sorted flag before doing this?
  #But samtools sort does not set it (not in the spec)!
  #samtools view -H unsort.bam
  #@HD    VN:1.0    SO:unsorted
  #samtools view -H sort.bam
  #@HD    VN:1.0    SO:coordinate
  #We could add it here, but VN is mandatory and we don't know the version of the sam format being used?
  #bwa doesn't seem to output the HD field, not do the docs suggest which spec is used for a given version
  #mailed Heng Lee regarding this
    
  #warn $cmd;  
  run_system_cmd($cmd);

#samtools view -h  /lustre/scratch110/ensembl/funcgen/alignments/homo_sapiens/GRCh37/ENCODE_Broad/HepG2_H2AZ_ENCODE_Broad.samse.bam -F 4 |  grep -vE '^[^[:space:]]+[[:blank:]][^[:space:]]+[[:blank:]][^[:space:]]+:[^[:space:]]+:MT:'  | grep -v '^MT' | grep -v '^chrM' | samtools view -uShb - | samtools sort - /lustre/scratch110/ensembl/funcgen/alignments/homo_sapiens/GRCh37/ENCODE_Broad/HepG2_H2AZ_ENCODE_Broad.samse.sorted_24137
#samtools rmdup -s /lustre/scratch110/ensembl/funcgen/alignments/homo_sapiens/GRCh37/ENCODE_Broad/HepG2_H2AZ_ENCODE_Broad.samse.bam.sorted_10600 - | samtools view -h - > /lustre/scratch110/ensembl/funcgen/alignments/homo_sapiens/GRCh37/ENCODE_Broad/HepG2_H2AZ_ENCODE_Broad.samse.filtered.sam
  #Add a remove duplicates step
  #-s single end reads or samse (default is paired, sampe)
  
  #todo, need to check infile and outfile aren't the same
  #else throw or output to tmp before mving back?
  
  $cmd  = "samtools rmdup -s ${sorted_bam_prefix}.bam - | ".   
    "samtools view -h${out_flag} - > $out_file";

  #Alternative with no rmdup...
  #$command .= $self->_bin_dir()."/samtools view -h ${input}_tmp.bam | gzip -c > $output";
  #warn $cmd;
  run_system_cmd($cmd);
  run_system_cmd("rm -f ${sorted_bam_prefix}.bam"); 
  
  return $out_file;
}


sub species_chr_num{
    my ($species, $val) = @_;

    ($species = lc($species)) =~ s/ /_/;

    my %species_chrs = (
                        homo_sapiens => {(
                                          'x' => 23,
                                          'y' => 24,
                                          'mt' => 25,
                                         )},

                        mus_musculus => {(
                                          'x'  => 20,
                                          'y'  => 21,
                                          'mt' => 22,
                                           )},

                        rattus_norvegicus =>  {(
                                                'x'  => 21,
                                                'y'  => 22,
                                                'mt' => 23,
                                               )},
                       );

    die("species not defined in chromosome hash") if(! exists $species_chrs{$species});

    return (exists $species_chrs{$species}{lc($val)}) ? $species_chrs{$species}{lc($val)} : $val;
}



#migrate this data to defs file!!??
#must contain all E! species and any other species which are used in local DB extractions
#NEED TO ADD FLY!!

sub species_name{
  my($species) = @_;
  my %species_names = (
		       "HOMO_SAPIENS", "human",
		       "MUS_MUSCULUS", "mouse",
		       "RATTUS_NORVEGICUS", "rat",
		       "CANIS_FAMILIARIS", "dog",
		       "PAN_TROGOLODYTES", "chimp",
		       "GALLUS_GALLUS", "chicken",
		       "SACCHAROMYCES_CEREVISIAE", "yeast",
		       "HUMAN",  "HOMO_SAPIENS",
		       "MOUSE", "MUS_MUSCULUS",
		       "RAT","RATTUS_NORVEGICUS",
		       "DOG", "CANIS_FAMILIARIS",
		       "CHIMP", "PAN_TROGOLODYTES",
		       "CHICKEN", "GALLUS_GALLUS",
		       "YEAST", "SACCHAROMYCES_CEREVISIAE",
		      );

  return $species_names{uc($species)};
}



#These subs are useful for implementing
#a farm mode in a run script, where a script can
#submit itself to the farm as slice based jobs

#strip cmd line params and associated arguments from a list
#should not be used to remove flag options i.e. no following args
#as this may cause removal of any following @ARGV;
#Can this be used on flattened args hash?

sub strip_param_args{
  my ($args, @strip_params) = @_;

  my $param_name;
  my $seen_opt = 0;

  foreach my $i(0..$#{$args}){

	if($args->[$i] =~ /^[-]+/){
	  $seen_opt = 0;#Reset seen opt if we seen a new one

	  ($param_name = $args->[$i]) =~ s/^[-]+//;

	  if(grep/^${param_name}$/, @strip_params){
		$seen_opt = 1;
	  }
	}

	#$args->[$i] = '' if $args->[$i] =~ /^[-]+farm/;#Only remove current flag
	#$seen_opt = 1 if $args->[$i] =~ /^[-]+skip_slices/;
	#$seen_opt = 1 if $args->[$i] =~ /^[-]+slice/;#Don't have full param name incase we have just specified -slice

	$args->[$i] = '' if $seen_opt;#Remove option and args following option
  }

  return $args;
}


sub strip_param_flags{
  my ($args, @strip_params) = @_;

  my @args = @$args;

  foreach my $flag(@strip_params){
	@args = grep(!/[-]+${flag}$/, @args);
  }

  return \@args;
}


sub url_from_DB_params {
  my $db_params = $_[0];
  
  my $host   = (exists $db_params->{'-host'})   ? $db_params->{'-host'} :
    throw('Missing mandatory -host paramter in DB parameters hash');
  my $user   = (exists $db_params->{'-user'})   ? $db_params->{'-user'} :
    throw('Missing mandatory -user paramter in DB parameters hash');
  my $dbname = (exists $db_params->{'-dbname'})   ? $db_params->{'-dbname'} :
    throw('Missing mandatory -dbname paramter in DB parameters hash');
  my $driver = (exists $db_params->{'-driver'}) ? $db_params->{'-driver'} : 'mysql';
  my $port   = (exists $db_params->{'-port'}) ? $db_params->{'-port'} : '';
  my $pass   = (exists $db_params->{'-pass'}) ? $db_params->{'-pass'} : '';
   
  return $driver.'://'.$user.':'.$pass.'@'.$host.':'.$port.'/'.$dbname;
}


sub validate_path{
  my ($path, $create, $is_dir, $alt_root) = @_;
            
  #how will alt_root and create interact
  #should they be mutually exclusive
  #or should we try both, but only create on the first root?
  #alt root will only work if path is a ref!
  #Is this the right way to do this, 
  #or should we handle root subbing in hive where we have access to both roots?
  
  #We should always test the dir first?
  #i.e. we don't want to mix files from two different dirs?
  
  #How will thi swork when we have set teh alignment dir
  #but we can't find something in there
  #we would also need to set alt alignment dir?
  #Whenever we set a dir, we shoudl also test and set an alt dir
  #if the first dir is absent, but we find the alt dir, then we should
  #set the alt as the main dir, then undef the alt dir
  
  #when ever we validate a path, we will then have access to an alt dir for 
  #which ever root we are using, or not if it doesn't exist
  
  #Still a little tricky about using files from two different dirs
  #we definately don't want to do this!
  #hard to manage this
  #maybe we just want to do this for fastqs and other 'reference' files etc?
  #all other intermediate piepline outptu/input should use the same root!
       
            
  if(! defined $path){
    throw('A path argument must be defined');
  } 
      
  #$label = ($label) ? ' '.$label : '';
  my $path_ref = ref($path);
  
  ### Build path  
  if($path_ref ne ''){
    
    if($path_ref ne 'ARRAY'){
      throw("Path argument is not valid:\t$path\n".
        'It must be a String or an Arrayref of string to pass to File::Spec'); 
    }
   
    #This works for files too instead of catfile
    #Will not prepend / if absent
    $path = File::Spec->catdir(@$path);
  }
  else{
    #Tidy up paths from config
    #This won't tidy up paths already stored in meta from init_pipeline.pl
    $path = File::Spec->canonpath($path); 
  } 
   
  ### Validate/create file/dir
  if(-e $path){
    
    if($is_dir){
      if( ! -d $path){
        throw("Found path but is not a directory:\n\t$path"); 
      }
      elsif($create && 
          ! -w $path){ #This chsould be writeable if we want to create it
        throw("Found path but is not writeable:\n\t$path");
      }
    }
    elsif(! -f $path){
      throw("Found path but is not a file:\n\t$path");
    }
    #could handle file create here?
  }
  elsif($is_dir && $create){
      #print STDOUT "Creating${label}:\n\t$path\n";
      
    eval { make_path($path) };
      
    if($@){
      throw("Unable to create:\n\t$path\n$@"); 
    }    
  }
  else{
    throw("Could not find:\n\t$path");
  } 
  
  return $path;  
}

#Funcgen added methods



1;
