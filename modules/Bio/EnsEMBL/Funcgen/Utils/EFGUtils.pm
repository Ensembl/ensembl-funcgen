=head1 LICENSE

Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


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

use warnings;
use strict;

use Digest::MD5;
use Bio::EnsEMBL::Utils::Exception qw( throw      );
use Bio::EnsEMBL::Utils::Scalar    qw( assert_ref );
use Scalar::Util                   qw( blessed    );
use File::Path                     qw( make_path  );
use File::Basename                 qw( dirname fileparse );
use File::Spec;
use Time::Local;
use FileHandle;
use Carp;

use base qw( Exporter );
use vars   qw( @EXPORT_OK );

@EXPORT_OK = qw(
  add_external_db
  add_hive_url_to_meta
  assert_ref_do
  assert_refs
  assert_refs_do
  backup_file
  convert_bam_to_sam
  create_Storable_clone
  file_suffix_parse
  generate_slices_from_names
  get_current_regulatory_input_names
  get_date
  get_files_by_formats
  get_file_format
  get_month_number
  get_my_scalar_name
  get_feature_file
  gunzip_file
  is_bed
  is_gzipped
  is_sam
  mean
  median
  open_file
  parse_DB_url
  path_to_namespace
  process_bam
  process_sam_bam
  run_backtick_cmd
  run_system_cmd
  scalars_to_objects
  set_attributes_by_my_scalar_names
  species_chr_num
  species_name
  strip_param_args
  strip_param_flags
  url_from_DB_params
  validate_checksum
  validate_path
  write_checksum
  );




#Split out methods into FileUtils.pm?


my %bed_strands = ('+' =>  1,
                   '-' => -1,
                   '.' =>  0);



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


=head2 create_Storable_clone

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


#$formats should be in preference order? Although this doesn't break things, it will just return a non-optimal file format

#Slightly horrible method to manage acquisition and conversion of files between
#supported formats (not necessarily feature formats)
#all_formats is necessary such that we don't redundant process files which are on the same conversion path when we have filter_format set


#There is a possibility that the formats provided might not have the same root, and so
#filter_from_format may be invlaid for one
#In this case two method calls might be require, hence we don't want to throw here if we can't find a file

  #This seems over-engineered! But we definitely need the ability to request two formats at the same time
  #to prevent parallel requests for the same file

  #Filtering will normally be done outside of this method, by the alignment pipeline
  #however, we must support it here incase we need to refilter, or we get alignment files
  #supplied outside of the pipeline


#what about if we only have the unfiltered file
#but we ask for filtered
#should we automatically filter?
#should we move handling 'unfiltered' to here from get_alignment_files_by_InputSet_formats?
#Is this too pipeline specific?
#what if some files don't use the 'unfiltered' convention?
#then we may get warnings or failures if the in_file and the out_file match
#would need to expose out_path as a parameter
#which would then need to be used as the in path for all subsequent conversion
#No this wouldn't work as it would change the in file to contain unfiltered
#which might not be the case.

sub get_files_by_formats {
  my ($path, $formats, $params) = @_;
  assert_ref($formats, 'ARRAY');
  $params ||= {};
  assert_ref($params, 'HASH');
  $params->{sort}     = 1 if ! defined $params->{sort};     #Always sort if we call process_$format
  #process_$format will never be called if $format file exists, hence no risk of a redundant sort
  #for safety, only set this default if filter_from_format is defined? in block below

  $params->{checksum} = 1 if ! defined $params->{checksum}; #validate and check

  if(scalar(@$formats) == 0){
    throw('Must pass an Arrayref of file formats/suffixes in preference order e.g. [\'sam\', \'bed\']');
  }

  my %conversion_paths = ( bam => ['bam'],
                           sam => ['bam', 'sam'],
                           bed => ['bam', 'sam', 'bed'],
                           #we always need the target format as the last element
                           #such that we can validate the filter_format e.g. for bam
                           #if the path array only has one element, it must match the key
                           #and this constitues calling filter_bam
                           #or if filter_format not set, just grabbing the bam file

                           #This approach prevents being able to | bam sort/filters through
                           #to other cmds, so may be slower if we don't need to keep intermediate files?

                           #Could also have non-bam rooted paths in here
                           #and maybe multiple path with different roots?
                         );

  my $can_convert           = 0;
  my $clean_filtered_format = 0;
  my $filter_format         = $params->{filter_from_format};
  my $all_formats           = $params->{all_formats};
  my $done_formats          = {};
  my ($files, @feature_files);

  #Add filter format if it is not in $formats
  if($filter_format &&
     (! grep(/^$filter_format$/, @$formats) ) ){
    unshift @$formats, $filter_format;
    $clean_filtered_format = 1;
  }

  #Attempt to get the first or all formats
  foreach my $format(@$formats){
    my $can_filter = 0;

    #Do this before simple file test as it is quicker
    if(grep(/^${format}$/, (keys %$done_formats))){ #We have already created this format
      next;
    }

    #Simple/quick file test first before we do any conversion nonsense
    #This also means we don't have to have any conversion config to get a file which
    
    #This is being undefd after we filter, so hence, might pick up a pre-exising file!
    if(! defined $filter_format){

       if(my $from_path = check_file($path.'.'.$format, 'gz', $params)){#we have found the required format
          $done_formats->{$format} = $from_path;
          next;
       }
    }


    ### Validate we can convert ###
    if(exists $conversion_paths{$format}){
      $can_convert = 1;

      if(defined $filter_format){

        if( ($conversion_paths{$format}->[0] ne $filter_format) &&
            ($all_formats) ){
          throw("Cannot filter $format from $filter_format for path:\n\t$path");
        }
        elsif((scalar(@{$conversion_paths{$format}}) == 1 ) ||
              (! $clean_filtered_format)){
          my $filter_method = 'process_'.$filter_format;
          $can_filter = 1;

          #Sanity check we can call this
          if(! ($filter_method = Bio::EnsEMBL::Funcgen::Utils::EFGUtils->can($filter_method))){
            throw("Cannot call $filter_method for path:\n\t$path\n".
              'Please add method or correct conversion path config hash');
          }

          #Set outfile here so we don't have to handle unfiltered in process_sam_bam
          #don't add it to $params as this will affected all convert methods
          (my $outpath = $path) =~ s/\.unfiltered$//o;

          #$format key is same as first element

          $done_formats->{$format} = $filter_method->($path.'.'.$filter_format, {%$params, 
                                                                          out_file => $outpath.'.'.$filter_format} );       
          #so we don't try and refilter when calling convert_${from_format}_${to_format}
 
          #delete $params->{filter_from_format};#Is this right?

          undef $filter_format; #Just for safety but not strictly needed
          $path = $outpath;

        }
      }
    }
    elsif($all_formats){
      throw("No conversion path defined for $format. Cannot acquire $format file for path:\n\t$path\n".
        'Please select a supported file format or add config and conversion support methods for $format');
    }

    ### Attempt conversion ###
    if($can_convert){
      #This now assumes that if $filter_format is set
      #convert_${filter_format}_${to_format} provides filter functionality

      if(scalar(@{$conversion_paths{$format}}) != 1){      #already handled process_${format} above
        #Go through the conversion path backwards
        #Start at last but one as we have already checked the last above i.e. the target format
        #or start at 0 if we have $filter_format defined
        my $start_i = (defined $filter_format) ? 0 : ($#{$conversion_paths{$format}} -1);

        for(my $i = $start_i; $i>=0; $i--){
          my $from_format = $conversion_paths{$format}->[$i];

          #Test for file here if we are not filtering! Else we will always go through
          #other formats and potentially redo conversion if we have tidied intermediate files
          if( (! defined $filter_format) &&
              (! grep (/^${from_format}$/, keys(%$done_formats)) ) ){
            my $from_path = $path.".${from_format}";

            if($from_path = check_file($from_path, 'gz', $params)){#we have found the required format
              $done_formats->{$from_format} = $from_path;
              #next; #next $x/$to_format as we don't want to force conversion
            }
          }


          #find the first one which has been done, or if none, assume the first is present
          if( (grep (/^${from_format}$/, keys(%$done_formats)) ) ||
              $i == 0){
            #then convert that to the next, and so on.
            for(my $x = $i; $x < $#{$conversion_paths{$format}}; $x++){
              my $to_format   = $conversion_paths{$format}->[$x+1];
              $from_format    = $conversion_paths{$format}->[$x];
              my $conv_method = 'convert_'.$from_format.'_to_'.$to_format;

              #Sanity check we can call this
              if(! ($conv_method = Bio::EnsEMBL::Funcgen::Utils::EFGUtils->can($conv_method))){
                throw("Cannot call $conv_method for path:\n\t$path\n".
                  'Please add method or correct conversion path config hash');
              }


              $done_formats->{$to_format} = $conv_method->($path.'.'.$from_format, $params);

              #Remove '.unfitlered' from path for subsequent conversion
              if(($i==0) &&
                defined $filter_format){
                $path =~ s/\.unfiltered$//o;
              }
            }

            last; #We have finished with this $conversion_path{format}
          }
        }


        if($clean_filtered_format && ($format eq $filter_format)){
          #filter_format is not our target format, so we need to keep going
          next; #$format
        }
        elsif(! $all_formats){  #else we have found the most preferable, yay!
          last;  #$format
        }
      }
    }
  } #end foreach my $format


  #Now clean $done_formats

  if($clean_filtered_format){
   #actually delete filtered file here?
    delete $done_formats->{$filter_format};
  }

  foreach my $format(keys %$done_formats){
    #doesn't matter about $all_formats here

    if(! grep(/^${format}$/, @$formats)){
      delete $done_formats->{$format};
    }
  }

  #test we have somethign to return!?
  #if( scalar(keys %$done_formats) == 0 ){
  #  throw('Failed to find any '.join(', ', @$formats)." files for path:\n\t$path");
  #}
  #don't do this as we may want to test for a filtered file, before attempting a filter
  #from a different path
  #This is caught in get_alignment_files_by_InputSet_formats

  return $done_formats;
}



#Is validate_checksum going to have problems as files are gunzipped
#Should validate checksum also handle .gz files i.e. check for entry without .gz, gunzip and validate?
#Maybe all checksums should be done on gunzipped files


#DAMMIT! Part of the filtering is currently done in SAM!!!!
#Need to fix this so we can drop sam file completely.

#There is a danger that a filter_format maybe specified but a pre_process_method
#never get called. This will have to be handled in the first convert_method in the path
#but we could put a method check in place?


#All pre_process_$format methods need to handle filter_from_format
#and should faciliatate filter and sort functions
#Can we merge this with process_sam_bam?
#and maintain this as a simple wrapper process_bam, which somply sets output format
#then we can also have process_sam as another wrapper method
#This would mean moving $params support to sort_and_filter_sam(process_bam)
#and also filter_from_format support and generate_checksum

#No this will make unflitered naming mandatory for process_sam!

#Calling pre_process assumes we want to at least convert, filter or just sort
#Otherwise we can simply just use the file
#Need to support sort flag. We might not want to sort if we already have a sorted bam
#Always sort when filtering?
#

sub process_bam{
  my ($bam_file, $params) = @_;
  $params ||= {};
  assert_ref($params, 'HASH');
  return process_sam_bam($bam_file, {%$params, output_format => 'bam'});
}

sub convert_bam_to_sam{
  my ($bam_file, $params) = @_;
  $params ||= {};
  assert_ref($params, 'HASH');
  return process_sam_bam($bam_file, {%$params, output_format => 'sam'});
}

#sub process_sam would need to check_file with gz suffix!


#Need to implement optional sort_and_filter_sam here?

sub convert_sam_to_bed{
  my ($sam_file, $params) = @_;
  my $in_file;

  if(! ($in_file = check_file($sam_file, 'gz', $params)) ){
    throw("Cannot find file:\n\t$sam_file(.gz)");
  }

  (my $bed_file = $in_file) =~ s/\.sam(\.gz)*?$/.bed/;
  run_system_cmd($ENV{EFG_SRC}."/scripts/miscellaneous/sam2bed.pl -uncompressed -1_based -files $in_file");

  if( (exists $params->{checksum}) && $params->{checksum}){
    write_checksum($bed_file, $params);
  }
  

  return $bed_file;
}


sub write_checksum{
  my ($file_path, $params) = @_;
  my ($signature_file, $digest_method);

  if(defined $params){
    assert_ref($params, 'HASH');
    $signature_file = $params->{signature_file};
    $digest_method  = $params->{digest_method};
  }

  $digest_method ||= 'hexdigest';
  my $md5_sig  = generate_checksum($file_path, $digest_method);
  my $file_name = fileparse($file_path);

  if(! defined $signature_file){
    $signature_file = $file_path.'.CHECKSUM';
  }

  my $checksum_row = $md5_sig."\t".$file_name."\t".$digest_method;

  #Update or create entry in signature_file
  my $sigfile_row;
  if(-f $signature_file){
    $sigfile_row = `grep '$file_name' $signature_file` ||
                    die("Failed to grep checksum signature file:\t$signature_file");
  }

  if($sigfile_row){ #Update entry
    my $sed_string = "sed -r 's/^.*[[:space:]]$file_name";

    #handle absent digest method here
    if($sigfile_row =~ /${digest_method}$/){
      $sed_string .= "[[:space:]]$digest_method";
    }

    run_system_cmd("${sed_string}\$/$checksum_row/' $signature_file > ${signature_file}.tmp; ".
                   "mv ${signature_file}.tmp $signature_file");
  }
  else{#create/append to file
    run_system_cmd("echo '$checksum_row' >> $signature_file");
  }

  return $md5_sig; #or return sig file?
}


sub generate_checksum{
  my ($file, $digest_method) = @_;
  
  if($file =~ /\.gz$/){
    throw("It is unsafe to generate checksums for a compressed file:\n\t$file");  
  }
  
  $digest_method ||= 'hexdigest';

  my $ctx = Digest::MD5->new;

  if(! $ctx->can($digest_method)){
    throw("Digest::MD5 cannot call method $digest_method, ".
     'please choose a valid digest method or omit for default hexdigest method');
  }

  open(FILE, $file) or throw("Cannot open file for md5 digest:\t$file");
  $ctx->addfile(*FILE);
  my $md5_sig = $ctx->$digest_method;
  close(FILE);

  return $md5_sig;
}

#assume the format of the file is:
#checksum_sig\tfilename\tdigestmethod


sub validate_checksum{
  my ($file_path, $params) = @_;
  my ($signature_file, $digest_method);

  if(defined $params){
    assert_ref($params, 'HASH');
    $signature_file = $params->{signature_file};
    $digest_method  = $params->{digest_method};
  }

  my $file_name = fileparse($file_path);

  if(! defined $signature_file){
    $signature_file = $file_path.'.CHECKSUM';
  }

  if(! -f $signature_file){
    throw("Failed to find checksum file:\t$signature_file\nPlease specify one as an argument, or create default file");
  }

  my $checksum_row = `grep -E '[[:space:]]$file_name\[[:space:]]*.*\$' $signature_file` ||
                      die("Cannot acquire $file_name checksum info from:\t$signature_file");
  my ($md5_sig, $sig_file_name, $sig_digest_method) = split(/\s+/, $checksum_row);


  if((! defined $sig_file_name) ||
     ($sig_file_name ne $file_name)){
    throw("Failed to find $file_name entry in checksum signature file:\n\t$signature_file");
  }

  #default to digest method in file
  $digest_method     ||= $sig_digest_method;
  #Also need to account for absent $sig_digest
  $sig_digest_method ||= $digest_method;


  if(defined $sig_digest_method){
    if($digest_method ne $sig_digest_method){
    throw("Specified digest method($digest_method) does not match method found in ".
      "checksum signature file($sig_digest_method):\n\t$file_path\n\t$signature_file");
    }
  }
  else{
    warn "Could not find digest method in checksum signature file, assuming $digest_method";
  }

  my $new_md5_sig = generate_checksum($file_path, $digest_method);

  if($md5_sig ne $new_md5_sig){
    #This could be due to mismatched digest methods
    throw("MD5 ($digest_method) checksum does not match signature file for:".
      "\n\tFile:\t\t$file_path\nSignature file:\t$signature_file");
  }

  return $md5_sig;
}

#todo purge_checksum from file, as we are deleting the file, or replacing it with a gz file


#could try and get the file suffix first to be a bit more directed here
#especially if the format list grows

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
    $was_gzipped = 1;
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


#Have individual format methods, as this is likely to hit performance
#would be nive to have a no thro mode, which would return the original value

sub convert_strand_from_bed {
  my $strand = $_[0];

  if(! defined $strand){
   throw('Strand argument is not defined');
  }

  if(! exists $bed_strands{$strand}){
    throw('Strand argument is not a valid bed strand, should be one of: '.
      join(' ', keys %bed_strands));
  }

  return $bed_strands{$strand};
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


=head2 path_to_namspace

  Arg [1]    : String - Module path
  Example    : my $path = '/my/module.pm';
               eval { require $path };

               if(! $@){
                 my $namespace = path_to_namespace($path);
                 my $obj       = $namespace->new();
               }

  Description: Simply converts a path to a namespace. Useful when you have used
               required a module dynamically rather than 'using' at compile time.
  Returntype : String
  Exceptions : Throws if arguments not defined
  Caller     : General
  Status     : At risk

=cut

sub path_to_namespace {
  my $path = $_[0];

  if(! defined $path){
   throw('No path argument defined');
  }

  $path =~ s/[\/]+/::/g;
  $path =~ s/\.pmc*$//;
  return $path;
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

 Example     : run_system_cmd($my_cmd);

 Exceptions  : Throws exception if system command returns none zero and no exit boolean not set

=cut

################################################################################

#Allow no_exit as some programs(tab2mage) give successful non-zero exit codes!
#Would be nice to catch warn output also, but system doesn not handle this well and would
#to redirect STDERR to a file to read back in. 

sub run_system_cmd{
  my ($command, $no_exit) = @_;
  my $exit_status = system($command);#don't nest below, as this may cause problems with $!
  return _handle_exit_status($command, $exit_status, $!, $no_exit);
}


################################################################################

=head2 run_backtick_cmd

 Description : Method to control the execution of the standard `backtick` command
               and return the results

 ReturnType  : Array or Scalar as required

 Example     : $Helper->debug(2,"dir=$dir file=$file");

 Exceptions  : throws exception if system command returns none zero

=cut

################################################################################

sub run_backtick_cmd{
  my $command = shift;
  
  my (@results, $result);
  
  if(wantarray){
    @results = `$command`;
  }
  else{
    $result  = `$command`;
  }
  
  _handle_exit_status($command, $?, $!); 
  return wantarray ? @results : $result;
}

sub _handle_exit_status{
  my ($cmd, $exit_status, $errno, $no_exit) = @_;
  
  my $exit_code = $exit_status >> 8; #get the true exit code
  
  if ($exit_status == -1) {
    warn "Failed to execute. Error:\t$errno\n";
  }
  elsif ($exit_status & 127) {
    warn sprintf("Child process died with signal %d, %s coredump\nError:\t$errno",
                 ($exit_status & 127),
                 ($exit_status & 128) ? 'with' : 'without');
  }
  elsif($exit_status != 0) {
    warn sprintf("Child process exited with value %d\nError:\t$errno\n", $exit_code);  
  }

  if ($exit_code != 0){

    if (! $no_exit){
      throw("System command failed:\t$cmd");
    }
    else{
      warn("System command returned non-zero exit code:\t$cmd");
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


#todo
# 1 add support for filter config i.e. which seq_regions to filter in/out
# 2 Sorted but unfiltered and unconverted files may cause name clash here
#   Handle this in caller outside of EFGUtils, by setting out_file appropriately
# 3 add a DESTROY method to remove any tmp sorted files which may persist after an
#   ungraceful exit. These can be added to a global $main::files_to_delete array
#   which should then also be undef'd in DESTROY so they don't persisnt to another instance

#This warning occurs when only filtering bam to bam:
#[bam_header_read] EOF marker is absent. The input is probably truncated.
#This is not fatal, and not caught. Does not occur when filtering with sort
#maybe we shoudl also be catchign $@ after samtools view -H $in_file?

sub process_sam_bam {
  my ($sam_bam_path, $params) = @_;
  my $in_file;

  if(! ($in_file = check_file($sam_bam_path, undef, $params)) ){
    throw("Cannot find file:\n\t$sam_bam_path");
  }

  $params ||= {};
  assert_ref($params, 'HASH');
  my $out_file      = $params->{out_file}           if exists $params->{out_file};
  my $sort          = $params->{sort}               if exists $params->{sort};
  my $filter_format = $params->{filter_from_format} if exists $params->{filter_from_format};
  my $fasta_fai     = $params->{ref_fai}            if exists $params->{ref_fai};
  my $out_format    = $params->{output_format}      if exists $params->{output_format};
  $out_format     ||= 'sam';
  my $fasta_fai_opt = '';

  if (defined $fasta_fai) {
    validate_path($fasta_fai);     $fasta_fai_opt = " -t $fasta_fai ";
    #-t spec is now optional, as not integrating the header into each file
    #is risky as it could get regenerated and be mismatched
    #which may produce unpredicatable behaviour or faults if the header
    #doesn't match the data.
  }

  #sam defaults
  my $in_format = 'sam';
  my $out_flag   = '';
  my $in_flag   = 'S';

  if($out_format eq 'bam'){
    $out_flag = 'b';
  }
  elsif($out_format ne 'sam'){
    throw("$out_format is not a valid samtools output format");
  }

  if($in_file =~ /\.bam$/o){     # bam (not gzipped!)
    $in_format = 'bam';
    $in_flag   = '';
  }
  elsif($in_file !~ /\.sam(\.gz)*?$/o){ # sam (maybe gzipped)
    throw("Unrecognised sam/bam file:\t".$in_file);
  }

  #This is odd, we really only need a flag here
  #but we already have the filter_from_format in the params
  if(defined $filter_format &&
     ($filter_format ne $in_format) ){
    throw("Input filter_from_format($filter_format) does not match input file:\n\t$in_file");
  }


  if(! $out_file){
    ($out_file = $in_file) =~ s/\.${in_format}(.gz)*?$/.${out_format}/;

    if(defined $filter_format){
      $out_file =~ s/\.unfiltered//o;  #This needs doing only if is not defined
    }
  }

  #Sanity checks
  (my $unzipped_source = $in_file) =~ s/\.gz$//o;
  (my $unzipped_target = $out_file)     =~ s/\.gz$//o;

  if($unzipped_source eq $unzipped_target){
    #This won't catch .gz difference
    #so we may have an filtered file which matches the in file except for .gz in the infile
    throw("Input and output (unzipped) files are not allowed to match:\n\t$in_file");
  }

  if($filter_format){

    if($in_file !~ /unfiltered/o){
      warn("Filter flag is set but input file name does not contain 'unfiltered':\n\t$in_file");
    }

    if($out_file =~ /unfiltered/o){
      throw("Filter flag is set but output files contains 'unfiltered':\n\t$in_file");
    }
  }
  elsif(! $sort &&
        ($in_format eq $out_format) ){
    throw("Parameters would result in no change for:\n\t$in_file");
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

  eval { $in_file_header =`samtools view -H $in_file` };

  #$! not $@ here which will be null string

  if($!){

    if(! defined $fasta_fai){
      throw("Could not find an in file header or a sam_ref_fai for:\t$in_file$!");
    }
  }
  elsif(defined $fasta_fai){
    #Now what do we do?
    #will -H now return the faig header, or the infile header?
    #I think the fai header (at least for bam)

    $fai_header = `samtools view -H $fasta_fai_opt $in_file`;

    if($!){
      throw("Failed to identify view valid header from:\t$fasta_fai\n$!");
    }

    #Just make sure they are the same
    if($fai_header ne $in_file_header){
      warn("CHANGE THIS TO A THROW WHEN WE HAVE FINISHED DEV! Found mismatched infile and fai headers for:\n\t$in_file\n\t$fasta_fai");
      #should this throw or just default to in file?
      #shoudl probably give the options to reheader and sub out this whole compare_sam_header thing
      #reheader mode, would that all old @SQ are present in new fai header?
      #can we get the header form the fai by simply providing an empty file?
      #cat '' | samtools view -HSt fasta.fai -

    }

    #Just undef this for simplicity
    $fasta_fai_opt = '';
  }


  my $cmd = "samtools view -h${in_flag} $fasta_fai_opt $in_file "; # Incorporate header into file

  if($filter_format){
    $cmd .= "-F 4 | ". #-F Skip alignments with bit set in flag (4 = unaligned)
    " grep -vE '^[^[:space:]]+[[:blank:]][^[:space:]]+[[:blank:]][^[:space:]]+\:[^[:space:]]+\:MT\:' ". #Filter MTs
    " | grep -v '^MT' | grep -v '^chrM' ";                                                              #Filter more MTs
  }

  #-u uncompressed bam (as we are piping)
  #-S SAM input
  #-t  header file (could omit this if it is integrated into the sam/bam?)
  #- (dash) here is placeholder for STDIN (as samtools doesn't fully support piping).
  #This is interpreted by bash but probably better to specify /dev/stdin?
  #for clarity and as some programs can treat is as STDOUT or anything else?
  #-b output is bam
  (my $tmp_bam = $in_file) =~ s/\.$in_format//;
  my $sorted_prefix = $tmp_bam.".sorted_$$";
  $tmp_bam .= ($sort) ? ".sorted_$$.bam" : '.tmp.bam';

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
  $cmd .= ' | samtools view -uShb - ';  #simply convert to bam using infile header
  $cmd .= ($sort) ? ' | samtools sort - '.$sorted_prefix : ' > '.$tmp_bam;
  #warn $cmd;
  run_system_cmd($cmd);

  #Add a remove duplicates step
  #-s single end reads or samse (default is paired, sampe)
  #Do this after alignment as we expect multiple reads if they map across several loci
  #but not necessarily at exactly the same loci which indicates PCR bias
  if($filter_format){
    $cmd = "samtools rmdup -s $tmp_bam - | ".
      "samtools view -h${out_flag} - > $out_file";
  }else{
    $cmd = "samtools view -h${out_flag} $tmp_bam > $out_file";
  }

  #warn $cmd;
  run_system_cmd($cmd);
  run_system_cmd("rm -f $tmp_bam");

  if( (exists $params->{checksum}) && $params->{checksum}){
    write_checksum($out_file, $params);
  }

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

#Very simple method to check for a file and it's compressed variants
#This may cause problems if the $file_path is already .gz



sub check_file{
  my ($file_path, $suffix, $params) = @_;

  my $found_path;

  if(-f $file_path){
    $found_path = $file_path;
  }
  elsif(defined $suffix){
    
    if($file_path =~ /\.${suffix}/o){
      $file_path =~ s/\.${suffix}$//;
      
      if(-f $file_path){
        $found_path = $file_path;  
      }
    }
    elsif(-f $file_path.'.'.$suffix){
      gunzip_file($file_path.'.'.$suffix);
      $found_path = $file_path;    
    }
  }
    

  if($found_path){
    my $validate_checksum;

    if(defined $params){
      assert_ref($params, 'HASH');
      $validate_checksum = (exists $params->{checksum}) ? $params->{checksum} : undef;
    }

    if($validate_checksum){
      validate_checksum($found_path, $params);
    }
  }

  return $found_path;
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

