# $Id: Collector.pm,v 1.7 2011-01-10 11:27:34 nj1 Exp $


=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

#Your Bio::Ensembl::Collection::Feature defs module should inherit from here
#This could be a local defs file which you have created and require'd into your script

#If your collections defs module refers to a Bio::EnsEMBL::Feature, 
#then it's adaptor should inherit from the collections defs module



package Bio::EnsEMBL::Funcgen::Collector;
#Move this to Bio::EnsEMBL::Utils::Collector for 59?


use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument  ('rearrange');
use Bio::EnsEMBL::Utils::Exception ('throw');
use Bio::EnsEMBL::Funcgen::ResultFeature;

#use base('Bio::EnsEMBL::Collection');#ISA

our ($pack_template, $packed_size, @window_sizes); #These get set in the FeatureAdaptor
#Make these constants and remove setter functionality in methods?
#Only really important for pack template and windows, maybe these if we are going to start


our $max_data_type_size = 16777216; #Default is 16MB for long blob
#we need to deduct the size of the rest of the record here!
#For a 2byte packet the smallest window size possible is:
#(slice->length/(16777216/2)
#so int(bin_size)+1
#Obviously have to use the largest slice here, for human chr1:
#249,250,621/(16777216/2) = 29.7???
#We may need to up this slightly to account for larger chrs?
#Implications on memory usage? Is it 4 times for blob manipulation?
#Does substr require this manipulation?
#This max_allowed_packet_size does not seem to translate directly to the size of the
#data being stored e.g. quite a bit more is needed.  ISG haven't got to the bottom of this yet.
#But have simply upped the config to 67108864 to handle the largest human chr.

our $max_view_width     = 500000;#Max width in Region In Detail;


#our %VALID_BINNING_METHODS
#Remove this in favour of can->('calculate_.$method) and coderefs?




#To do
# 1 DONE Merge in Collection code, (no need to do this, removed inheritance)
# 2 Write simple BED input to flat file output.
# 3 Separate store method so we can simply get, then wrap store around this
# 4 Test get method with slice adjusts
# 5 separate set_config? 
# 6 optimise generate_bin_chunks to handle just one window size for display?
# 7 Handle packed_size pack_template as methods  constants
# 8 Provide override method in basefeature adaptor which will use package constant in feature adaptor
# This is because these are really adaptor config, the collector only needs to know the
# packed_size, and in the absence of an feature adaptor also provides the default methods for both.
# If we substr in the API then we need to set sensible limits on blob size, otherwise we will have to unpack a lot of data
# to get at the slice we want.
# OR
# Change adaptor to substr in DB based on known blob ranges/window size
# and stitch together any which cross boundaries. This depends on speed of substr at end of large blob TEST!
# Load with current code first and test this before making either change!
# Delete empty (non-0) collections? i.e. For seq_regions which do not have any source features.
#
# 9 Handle PAR/HAP regions using fetch_normalised_slice_projections This has to be done in the feature adaptor! Then restrict  to non_dup regions in calling script



=head2 new

  Args       : None
  Example    : my $collector = Bio::EnsEMBL::(Funcgen|Compara|Variation::)Collector::FEATURE->new;
               $collector->store_windows_by_Slice($slice);
  Description: Simple new method to enable use of collector when not inherited by
               a descendant of Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor
  Returntype : Bio::EnsEMBL::Funcgen::Collector
  Exceptions : None
  Caller     : Collector script
  Status     : At Risk

=cut

sub new{
  return bless {}, $_[0];#Simple blesses this class as an empty hash

  #Do not set anything here
  #As will not be first in ISA for feature adaptors
  #Hence not guaranteed to be called

}

#Setter/Getter methods if we don't have a dedicated Collector module
#to set the package variables in? Also to allow overriding of defaults.
#This can be used by the write_collection method
#to determine when to build and store a compressed collection
#Effectively the max size of the data type you are using to store 
#a compressed score. defaults to max for long blob 16MB



#Generic method, but only ever called by write_collection in descendant

sub new_assembly{
  my ($self, $new_ass) = @_;

  if($new_ass){
	#Validate new assm to project to


	$self->{'new_assembly'} = $new_ass;
  }

  return $self->{'new_assembly'};
}


sub max_data_type_size{
  my ($self, $size) = @_;

  #Validate is sensible integer?

  if($size && ! int($size)){
	throw("max_data_type_size must be a integer of bytes, not $size");
  }
  elsif($size){
	$self->{'max_data_type_size'} = $size;
  }
  elsif(! defined $self->{'max_data_type_size'}){
	#default set at head of this module or in descendant Collector
	$self->{'max_data_type_size'} = $Bio::EnsEMBL::Funcgen::Collector::max_data_type_size;
  }

  return $self->{'max_data_type_size'};
}

sub max_view_width{
  my ($self, $size) = @_;

  #Validate is sensible integer?

  if($size && ! int($size)){
	throw("max_view_width must be a integer, not $size");
  }
  elsif($size){
	$self->{'max_view_width'} = $size;
  }
  elsif(! defined $self->{'max_view_width'}){
	#default set at head of this module or in descendant Collector
	$self->{'max_view_width'} = $Bio::EnsEMBL::Funcgen::Collector::max_view_width;
  }

  return $self->{'max_view_width'};
}


sub bins_per_record(){
#$collector_class::bins_per_record = ($collector_class::max_data_type_size/$collector_class::packed_size);#This should be done dynamically as we may redefine either of these variables?

  my ($self) = shift;

  return int($self->max_data_type_size/$self->packed_size);
}


#The defaults for these should be defined in the feature/format specific Collector descendant
#either by specifying the package variables or using config attrs to set methods?
#general config should be parsed here.
#rename bin_method?

sub bin_method{
  my ($self, $method) = @_;

  if($method || ! $self->{'bin_method'}){

	if($method){
	  $self->{'bin_method'} = $method;
	  #should test can here? or validate versus hash?
	}
	elsif(! $self->{'bin_method'}){
	  
	  if (! defined  $Bio::EnsEMBL::Funcgen::Collector::bin_method){
		throw('Must pass a bin_method in the config or define $Bio::EnsEMBL::Funcgen::Collector::bin_method in your Collector');
	  }
	  
	  $self->{'bin_method'} = $Bio::EnsEMBL::Funcgen::Collector::bin_method;
	}

	#or current validate method if we are keeping the method in the if/else block
	

  	#if(! $self->can("calculate_${method}"))){
	#throw("$method is no a valid a valid binning method");
	#}
  
  }

  return $self->{'bin_method'};
}

#We could replace this with a hash of bin_methods and models?
#This could then be used to validate
#Altho if we are going to commodotise the bin methods, then we need to be able to 
#define this in the child Collector. Could still do this by modifying the method/model
#hash from the child Collector

sub bin_model{
  my ($self, $bin_model) = @_;
  
  if($bin_model || ! $self->{'bin_model'}){

	if($bin_model){
	  $self->{'bin_model'} = $bin_model;
	}
	elsif(! $self->{'bin_model'}){
	  
	  #Set as global constant defined in descendant Collector
	  if (! defined  $Bio::EnsEMBL::Funcgen::Collector::bin_model){
		throw('Must pass -bin_model in the config or define $Bio::EnsEMBL::Funcgen::Collector::bin_model in your Collector');
	  }
	  
	  $self->{'bin_model'} = $Bio::EnsEMBL::Funcgen::Collector::bin_model;
	}
	
	#Need to validate bin models here
	throw('Bio::EnsEMBL::Funcgen::Collector does not yet support non-SIMPLE bin models')	if $self->{'bin_model'} ne 'SIMPLE';
  }

  return $self->{'bin_model'};
}

#This can be overridden by adaptor method
#At present this could cause problems as we can pass window sizes in the config, but they will never be set
#as adaptor method is not a setter. Adaptor method should throw if we try and set them as this could cause problems when fetching and not knowing the custom sizes?

sub window_sizes{
  my ($self, $sizes) = @_;

  if($sizes || ! $self->{'window_sizes'}){

	if($sizes){
	  $self->{'window_sizes'} = $sizes;
	}
	else{#! $self->{'windows_sizes'

	  if (! @window_sizes){
		throw('Must pass -windows_sizes in the config or define @Bio::EnsEMBL::Funcgen::Collector::window_sizes in your Collector');
	  }
	  
 	  @{$self->{'window_sizes'}} = @window_sizes;
	}
	
	if(ref($self->{'window_sizes'}) ne 'ARRAY' ||
	   scalar(@{$self->{'window_sizes'}}) == 0){
	  throw('window_sizes must be an arrayref of at least one window size');
	}
  }

  return $self->{'window_sizes'};
}




#Optional attrs dependant on whether Collection is packed
#Can be redefined in the adaptor but becareful never to redefine the actual values
#As these should really be constants for a given Collector
#What is best here? We only want pack methods for storing/fetching compressed collections
#Move this to base feature adaptor and define attrs as constants using
#package variable? Or directly in new?
#Then direct modification will be caught.
#Just leave here for now.

#Caller _obj_from_sth/store

sub pack_template{
  my ($self, $template) = @_;

  if($template){
	$self->{'pack_template'} = $template;
  }
  elsif(! $self->{'pack_template'}){

	#Set as global constant defined in descendant Collector
	
	if (! defined  $Bio::EnsEMBL::Funcgen::Collector::pack_template){
	  throw('Must pass a per score pack_template in the config or define $Bio::EnsEMBL::Funcgen::Collector::pack_template in your Collector');
	}

	$self->{'pack_template'} = $Bio::EnsEMBL::Funcgen::Collector::pack_template;
  }

  return $self->{'pack_template'};

}

#Caller _obj_from_sth/store & current_packed_size

sub packed_size{
  my ($self, $size) = @_;

  if($size){

	if(! int($size)){
	  throw("$size is not an integer, must pass a size integer for packed_size which specifies size of pack_template:\t".$self->pack_template);
	}

	$self->{'packed_size'} = $size;
  }
  elsif(! $self->{'packed_size'}){

	#Set as global constant defined in descendant Collector
	
	if (! defined  $Bio::EnsEMBL::Funcgen::Collector::packed_size){
	  throw('Must pass a packed_size(wrt to pack_template) in the config or define $Bio::EnsEMBL::Funcgen::Collector::packed_size in your Collector');
	}

	$self->{'packed_size'} = $Bio::EnsEMBL::Funcgen::Collector::packed_size;
  }

  return $self->{'packed_size'};

}

#These methods are used by the descendant Collector
#For caching infor whilst building collections
#This is used to log how big a collection has grown before storing

sub current_packed_size{
  my ($self, $wsize) = @_;
  
  #$self->{'current_packed_size'}{$wsize} ||= 0;

  #if(defined $cps){
#	$self->{'current_packed_size'}{$wsize} = $cps;
#  }
#  else{
#	return $self->{'current_packed_size'}{$wsize};
#  }

  return (scalar(@{$self->score_cache($wsize)})*$self->packed_size);

}


sub score_cache{
  my ($self, $wsize, $scores) = @_;

  $self->{'score_cache'}{$wsize} ||= [];

  if(defined $scores){
	push @{$self->{'score_cache'}{$wsize}}, @{$scores};
  }
  else{
	#Do this here to stop passing the ref everytime
	#Will this be faster?
	#Would certainly be faster if we were not returning a ref
	return $self->{'score_cache'}{$wsize};
  }
}

#These last methods are only used for the 0 wsize
#natural resolution and ar wrt the orig_slice passed
#to store_windows_by_Slice

sub collection_start{
  my ($self, $wsize, $sr_start) = @_;

  if(defined $sr_start){
	$self->{'collection_start'}{$wsize} = $sr_start;
  }
  else{
	return $self->{'collection_start'}{$wsize};
  }
}

sub collection_end{
  my ($self, $wsize, $sr_end) = @_;

  if(defined $sr_end){
	$self->{'collection_end'}{$wsize} = $sr_end;
  }
  else{
	return $self->{'collection_end'}{$wsize};
  }
}

sub collection_strand{
   my ($self, $wsize, $strand) = @_;

   if(defined $strand){
	 $self->{'collection_strand'}{$wsize} = $strand;
   }
   else{
	 return $self->{'collection_strand'}{$wsize};
  }
}



=pod

sub _create_feature {
  my ( $this, $feature_type, $args ) = @_;

  my $feature = $this->SUPER::_create_feature( $feature_type, $args );

  if ( !$this->_lightweight() ) {
    my ( $phase, $end_phase, $stable_id, $version, $created_date,
         $modified_date, $is_current )
      = rearrange( [ 'PHASE',        'END_PHASE',
                     'STABLE_ID',    'VERSION',
                     'CREATED_DATE', 'MODIFIED_DATE',
                     'IS_CURRENT'
                   ],
                   %{$args} );

    push( @{$feature},
          $phase, $end_phase, $stable_id, $version, $created_date,
          $modified_date, $is_current );
  }

  return $feature;
}

sub _create_feature_fast {
  my ( $this, $feature_type, $args ) = @_;

  my $feature =
    $this->SUPER::_create_feature_fast( $feature_type, $args );

  return $feature;
}


#This might not be sensible for Features which are split across tables

sub _tables {
  my ($this) = @_;

  my @tables = $this->SUPER::_tables();

  if ( $this->_lightweight() ) {
    return ( $tables[0] );
  }

  return @tables;
}


sub _columns {
  my ($this) = @_;

  my @columns = $this->SUPER::_columns();

  if ( $this->_lightweight() ) {

	#What is this doing?
	#Probably not sensible for ResultFeature
    @columns[ 5 .. $#columns ] = map( 1, 5 .. $#columns );
  }

  return @columns;
}

#Also not sensible for objects spread across several tables

sub _default_where_clause {
  my ($this) = @_;

  if ( $this->_lightweight() ) {
    return '';
  }

  return $this->SUPER::_default_where_clause();
}

=cut



#This need to be generic
#Again we need to pass an accessor method/reference?
#Will be some sort of generic fetch for feature adaptors
#or while loop for in flat file accessor
#rollback to be handled in caller?

# To do
#
# 1 Allow variable chunks lengths (so we only have one resolution of windows?)
#   This will allow SNP collections which currently define classification i.e colour
#   Density of SNPs within window will define shading. Count will be displayed in zmenu
#   This maybe something we have to do in the descendant
#
# 2 Implement collection param definition in/from descendant


# return collection config from adaptor fetch
# window size
# fixed width?
# render/collection style?
# This chould be implemented in BaseFeatureAdaptor::generic_fetch?
# Or could be done in the calling fetchmethod?
#

#need to change this to get_window_bin_by_Slice
#to enable generating bins on uncompressed data
#Need to remove all counts and store based code to store caller
#this would mean removing any pack based code too
#separate set_config method


#Probelm here is size of slice?
#We need to generate bins all in one go, but also need to store at interval
#so as not to explode memory
#Do we need to separate the window generation from the bin generation code?


#Define the optimal way to generate windowed data by
#finding the most common denominator

sub _define_window_chunks{
  my ($self, $window_sizes, $max_view_size) = @_;

  ### DEFINE CHUNKS WRT WINDOWS

  #Shortcut for on the fly uncompressed collection retrieval
  #if(scalar(@wsizes) = 1){
  #
  #}
  #else{
   
  #Calulate sensible slice length based on window sizes
  my @wsizes = sort {$a <=> $b} @$window_sizes;
  
  #We need a default when only calculating 0 resolution
  #Will binning code work with only 0 resolution?
  if((scalar(@wsizes) == 1) &&
	 $wsizes[0] == 0){
	return { $self->max_view_width => [0] };
  }
  

  my $multiplier = int($max_view_size/$wsizes[$#wsizes]);
  my $chunk_length = $multiplier * $wsizes[$#wsizes];
  my $not_divisible = 1;
  my %chunk_windows;#Registry of chunk lengths to run with windows
  my %workable_chunks = map {$_ => {}} @wsizes;
  delete $workable_chunks{'0'};#get rid of natural resolution as this will always work


  while($not_divisible && $chunk_length != 0){
	$not_divisible = 0;
	
	foreach my $wsize(@wsizes){
	  next if $wsize == 0;#Special wsize for normal data
	  
	  #Set not divisible if modulus is true
	  if($chunk_length % $wsize){
		$not_divisible = 1;
	  }
	  else{
		#No need to listref here?
		$workable_chunks{$wsize}{$chunk_length} = [];
	  }
	  
	  #warn "chunk length is $chunk_length and not_divisible is $not_divisible";
	}
	
	#Gradually shrink the length until we find a workable slice length for all windows
	$chunk_length -= $wsizes[$#wsizes] if $not_divisible;
  }
  
  my %chunk_sets;
  
  
  if($chunk_length == 0){
	print "Could not find chunk length for all window sizes, attempting to subset windows using alternate slice length\n";
	
	foreach my $wsize(keys %workable_chunks){
	  #Loop through windows, seeing if they are workable in the other windows
	  
	  foreach my $chunk(keys %{$workable_chunks{$wsize}}){
		
		foreach my $other_wsize(keys %workable_chunks){
		  next if $wsize == $other_wsize;
		  
		  if(exists $workable_chunks{$other_wsize}{$chunk}){
			#only push it onto the other wsize, as we will do the reverse later			
			$chunk_sets{$chunk}{$wsize} = undef;
		  }
		}
	  }
	}
	
	#Now we have a register of co-occurence of wsizes with repect to chunks
	#Loop through finding the least amount of sets with the longest chunk length?
	#There is no way to decide which is best?
	#we could calculate the number of loops? Factored by the chunk length?
	#Let's just print out and see what we get
	
	#warn "chunk sets are :\n".Data::Dumper::Dumper(\%chunk_sets);


	#For now let's just take the one which has the most windows and the longest chunk
	#Then we just get the largest which handles the rest.

	#define possible set lengths
	my $i = 0;
	my %set_lengths;
	map {$set_lengths{$i} = []; $i++} @wsizes;
	delete $set_lengths{'0'};#get rid of natural resolution as this will always work


	#store chunks lengths for each set size
	foreach my $chunk(keys %chunk_sets){
	  my $set_size = scalar(values %{$chunk_sets{$chunk}});
	  push @{$set_lengths{$set_size}}, $chunk;
	}		
  
	#Now we get the biggest set with the longest length;
	my $largest_size = scalar(@wsizes);#scalar here as we are disregarding natural resolution of 0 in loop
	my $found_largest_set = 0;

	while(! $found_largest_set){
	  $largest_size--;

	  if(scalar(@{$set_lengths{$largest_size}}>0)){
		$found_largest_set = 1;
	  }
	}


	#We should be able to loop this bit, to find all the biggest sets.
	my ($largest_chunk) = sort {$b<=>$a} @{$set_lengths{$largest_size}};
	#we could even be selective here, but let's just take the first one for now

	my @largest_windows = keys %{$chunk_sets{$largest_chunk}};
	@{$chunk_windows{$largest_chunk}} = @largest_windows;
	print "Largest chunk $largest_chunk($largest_size) contains windows: @largest_windows\n";

	my %remaining_windows = map {$_ => {}} @wsizes;
	delete $remaining_windows{'0'};#get rid of natural resolution as this will always work
	map { delete $remaining_windows{$_} } @largest_windows;	
	my $remaining_set_size = scalar(keys %remaining_windows);

	#swapping to array here for practicality, would need to maintain hash if we need to iterate
	my @rwindows = keys %remaining_windows;

	#This could just be one window, but this will not be inthe co-occurence hash %chunk_sets
	#Hence the normal approach will not work. and we just want to find a suitably large chunk for this one window.
	my $next_chunk;

	if(scalar(@rwindows) == 1){
	  #we just want to find a suitably large chunk for this one window.
	  my ($last_window) = @rwindows;

	  $multiplier = int(500000/$last_window);
	  $next_chunk = $multiplier * $last_window;
	}
	else{
	  #Now were are doing something very similar to above
	  #populating a set_size chunk length registry
	  #my %seen_hash;
	 	  

	  foreach my $chunk(sort {$b<=>$a} @{$set_lengths{$remaining_set_size}}){
	  	my $seen_count = 0;
		
		foreach my $rwindow(@rwindows){
		  
		  $seen_count++ if grep/$rwindow/, (values %{$chunk_sets{$chunk}});
		}
		
		if ($seen_count == $remaining_set_size){
		  
		  $next_chunk = $chunk;
		  last;
		}
	  }
	}

	@{$chunk_windows{$next_chunk}} = @rwindows;


	if($next_chunk){

	  print "Found next chunk length $next_chunk contains remaining windows:\t@rwindows\n";

	#Now we want to cycle through all the set lengths which could contain the ones not in the first
	#so we need to
	}
	else{
	  warn "Need to write iterative sub for set definition";
	  throw('Could not find workable slice length for remaining windows: '.
			join(', ', @rwindows));

	}
  }
  else{
	@{$chunk_windows{$chunk_length}} = keys(%workable_chunks);
   	print "Found workable chunk length($chunk_length) for all window sizes:\t".
	  join(' ', @{$chunk_windows{$chunk_length}})."\n";
  }

  return \%chunk_windows;
}


#Let's concentrate on store function first before we split out into store and fetch methods
#How will this work with the Bed parser?
#The descendant collector will sort the input and detect the current slice before calling 
#store_window_bins_by_Slice.  This may require some caching of line or seeking as we will see the next slice before we have a chance to set it.
#This will store as ResultFeature collections, so maybe we need to separate the input from output code?
#i.e. Bed parser/wrapper
#     ResultFeatureAdaptor wrapper
#These 

#Problem with passing window_sizes here
#We need to check that they aren't already defined a class variables as this could potentially
#screw up retrieval, expect for only 0 or all but 0
#Should we remove this config and force the class variable to be set in the 'adaptor'
#Method is then only used internally, make private or only getter? Set by changing class vars?

sub store_window_bins_by_Slice{
  my ($self, $slice, %config) = @_;

  my ($window_sizes, $logic_name, $bin_method, $fetch_method_ref, $max_view_width, 
	  $max_data_type_size, $pack_template, $packed_size, $bin_model, $new_assm, $skip_zero_window) =
    rearrange( [ 'WINDOW_SIZES', 'LOGIC_NAME', 'BIN_METHOD', 'FETCH_METHOD_REF', 'MAX_VIEW_WIDTH', 'MAX_DATA_TYPE_SIZE', 'PACK_TEMPLATE', 'PACKED_SIZE', 'BIN_MODEL', 'NEW_ASSEMBLY', 'SKIP_ZERO_WINDOW'], %config );

  warn "Need to be careful here about cleaning start end strand caches between serially run slices";



  ### VAILDATE VARS/CONFIG

  #This could be done once in set_config, could then remove setter bahviour from attr methods?
  #All default defs params/methods can be overridden by config params
  #Attrs used in this method
  $bin_method   = $self->bin_method($bin_method);
  $bin_model    = $self->bin_model($bin_model);
  #$window_sizes = $self->window_sizes($window_sizes);#Now done below
  #Set to undef if we ave empty array
  $window_sizes = undef if (ref($window_sizes) eq 'ARRAY' && scalar(@$window_sizes) == 0);
  #Attrs used in other (store) methods
  $self->pack_template($pack_template);
  $self->packed_size($packed_size);
  $self->max_data_type_size($max_data_type_size);
  $self->max_view_width($max_view_width);

  #Other vars
  $self->new_assembly($new_assm);

  #Need to validate slice here

  warn "temp hack for bin_method validation";
  $bin_method = $self->validate_bin_method($bin_method);
    
  ### Set window_sizes
  
  if($self->new_assembly){
	print "Assembly projection may cause problems for large Collections, defaulting to window_sizes = (0)\n";
	#Then build the bins on the projected 0 level single ResultFeatures

	#Test we haven't explicity set window_sizes to be soemthing else

	if($window_sizes &&
	   ! ( scalar(@$window_sizes) == 1 && $window_sizes[0] == 0)){

	  throw("You have set window_sizes config which are not safe when projecting to a new assembly($new_assm), please omit window_sizes config or set to 0");

	}
	$window_sizes = $self->window_sizes([0]);
  }
  else{

	if($window_sizes && $skip_zero_window && grep/^0$/,@$window_sizes){
	  throw("You have specied skip_zero_window and window_size 0 in your config, please remove one of these");
	}
	elsif($window_sizes && ! grep/^0$/,@$window_sizes){
	  $skip_zero_window = 1;
	  unshift @$window_sizes, 0;#re-add 0 window as we need this to build the collections
	}

	$window_sizes = $self->window_sizes($window_sizes);
  }
  

  #This is already done in the script
  if($skip_zero_window && $new_assm){
	throw("You cannot -skip_zero_window or omit 0 from -window_sizes when projecting to a new assembly($new_assm) which should only be generated using window_size=0");
  }
 


  ### Rollback previously stored features
  
  if($self->can('rollback_Features_by_Slice')){
	$self->rollback_Features_by_Slice($slice);
  }
  else{
	#This is currently the only warn output we can't get rid off
	warn ref($self)." cannot rollback_Features_by_Slice. This may result in duplicate Collections being stored if there is pre-existing data";
  }



  ### PROCESS CHUNKS

  #Not lightweight as we will be storing them
  # Temporarily set the collection to be lightweight???
  #my $old_value = $this->_lightweight();
  #if   ( defined($lightweight) ) { $this->_lightweight($lightweight) }
  #else                           { $this->_lightweight(1) }

  my %chunk_windows = %{$self->_define_window_chunks($self->window_sizes, $self->max_view_width)};
  my (%counts, $store_natural);
  $store_natural = grep/^0/, @$window_sizes;
  $counts{0}=0;#Set natural res count to 0
  my $slice_end        = $slice->end;
  my $orig_slice       = $slice;
  my $orig_start       = $slice->start;
  #my $slice_adj = $slice->start - 1;#Removed this as we are now generating features local to orig_slice
  #start/end conversion will be done in write/store_collection
  my $region           = $slice->coord_system_name;
  my $version          = $slice->coord_system->version;
  my $seq_region_name  = $slice->seq_region_name;
  my $strand           = $slice->strand;
  my $only_natural     = 0;
  #my $slice_adj = 0;


  #We need to account for only 0 here when doing projection
  #The chunk window is set to max_view_widht in _define_chunk_windows

  $only_natural = 1 if $store_natural && scalar(@$window_sizes) == 1;
  $store_natural = 0 if $skip_zero_window;
  #SHould really test these two, but should already be caught by now

  #Set the initial  collection_start to orig_start
  #Could default to 1, but we may not be starting from 1
  #This is not the case for 0 wsize where it must always be 
  #The first feature start


  for my $wsize(@{$self->window_sizes}){

	next if $wsize == 0;# && $skip_zero_window;#We never want to assume start of 0 window collection
	$self->collection_start($wsize, $orig_start);
  }


 
  foreach my $chunk_length(sort keys %chunk_windows){

	print "Processing windows ".join(', ', @{$chunk_windows{$chunk_length}}).
	  " with chunk length $chunk_length\n";
	map  $counts{$_} = 0, @{$chunk_windows{$chunk_length}};	#Set window counts to 0

	#Now walk through slice using slice length chunks and build all windows in each chunk
	my $in_slice     = 1;
	my $start_adj    = 0;
	my ($sub_end, $features, $bins);
	my $sub_start    = 1;
	my $slice_length = $slice->length;
	

	#Can we subslice and then exclusivly use bin_start(local to orig_slice)
	#Then we never have to deal with sr coord until we store
	#This should me we never have to do the sr conversion unless we 
	#use a slice which doesn't start at 1(PAR or test)
	#Always create in local coords for fetch
	#Then change to seq_region coords for store if required


	while($in_slice){
	  #$sr_start       = $slice_start + $start_adj;
	  $sub_start += $start_adj;

	  #$slice_start = $sr_start;#Keep for next slice
	  #$sr_end   = $sr_start + $chunk_length - 1;
	  $sub_end   = $sub_start + $chunk_length - 1;

	  #Last chunk might not be the correct window length
	  #Hence why we should do this on whole chromosomes
	  if($sub_end >= $slice_length){
		#$sub_end = $slice_end;
		#No longer set to slice end, as we don't want to corrupt the bin definition?
		#Surplus bins are removed in store/write_collection in caller
		#We could simply add the largest window the the end of the slice?
		#Then we will only build the minimum of excess bins?
		#This should be okay for bin calcs
		#But may screw up bin trimming in caller as we currently expect $ub_end to be a valid bin end
		#for all wsizes
		#bin trimming should handle this, but this will corrupt the bin definition???
		#bin definition is depedant on method
		#So this method need to be agnostic
		#And deal with the rest in descendant
		$in_slice = 0;
	  }
	  
	
	  $slice = $slice->adaptor->fetch_by_region($region, $seq_region_name, ($sub_start + $orig_start -1), ($sub_end + $orig_start - 1), $strand, $version);
	  #Can't subslice as this will not clip if we go over the length of the slice, unlike normal slice fetching
	  #hence we cannot rely on this
	  #$slice = $orig_slice->sub_Slice($sub_start, $sub_end, $orig_slice->strand);
	  #warn "got sub slice $slice as $sub_start - $sub_end from ".$orig_slice->name;

	  
	  ### Grab features and shift chunk coords
	  #features may already be a 0 wsize collection if we have projected from an old assembly
	  #Could move this check to get_Features_by_Slice?
	  
	  #e.g. [ $features, \%config ]
	  $features = $self->get_Features_by_Slice($slice);
	  #next if scalar(@$features) == 0;#We want to store values for all windows

	  if( (@$features) &&
		  (ref($features->[0]) =~ /Bio::EnsEMBL::Funcgen::Collection/) ){#Change to isa 'Bio::EnsEMBL::Collection
		
		#Check that the returned feature/collections support window_size

		if($features->[0]->can('window_size')){
		  
		  if($features->[0]->window_size != 0){
			throw("You are trying to generated Collections from a non-zero window sized Collection:\t".$features->[1]->{'window_size'});
		  }
		  
		  #This should never happen
		  if(! $skip_zero_window){ 
			throw('You have retrieved data from a Collection which without using -skip_zero_window i.e. you are trying to generate overwrite the data you are generating the Collections from');
		  }
		}
		else{
		  throw('Something si wrong, the Collection you have retrieved does not support the method window_size');
		}
	  }



	  #Set collection start here for 0 window_size
	  if(@$features && $store_natural && ! defined $self->collection_start(0)){
		$self->collection_start(0, ($features->[0]->start + $sub_start));
	  }

	  

	  $start_adj = $chunk_length if($in_slice);
	  	  
	
	  
	  #This should return a hash of window size => bin array pairs
	  if(! $only_natural){
		
		$bins = $self->_bin_features_by_window_sizes(
													 -slice         => $slice,
													 -window_sizes  => $chunk_windows{$chunk_length},
													 -bin_method        => $bin_method,
													 -features      => $features,
													);


		

	  }

	  #my $bin_start = $sr_start + $slice_adj;#This was only required for storing individual bins
	  #Could calc bin_start + slice_adjust ahere for all features
	  #Doing this will break old code for single window collections
	  
	  #This is sr start and should be local to orig_slice!
	  



	  #We need to handle strandedness of slice!?	  
	  
	  #Store all normal features in result_feature


	  
	  if($store_natural){
		
		foreach my $feature(@$features){
		  $counts{0}++;
		  #warn "storing ".join(', ', ($feature->start, $feature->end, $feature->strand, $feature->scores->[0]));
		  

		  #Should we handle bin trimming here for overhanging slices
		  #Then counts wil be correct and wont have to do in caller

		  #We could stop here if the feature seq_region start > orig_slice end
		  #Current done in write/store_collection
		  #This may mean working in seq_region values rather than slice values


		  #write_collection is implemented in descendant e.g. Bio::EnsEMBL::Funcgen::Collector::ResultFeature
		  #as wrapper to adaptor store method or print to file
		  
		  #These params need to be generated in a way defined by the descendant
		  #
		  
		  if($bin_model eq 'SIMPLE'){
			#We need to pass the slice with this so we can sub slice when storing
			#the collection and set the start/end to 1 and length of slice
			#we still need to store the first start to be able to sub slice correctly
			
			$self->collection_start(0, ($feature->start + $sub_start));

			#Need to pass strand for 0 resolution
			$self->write_collection(0,
									$orig_slice,
									#These are now wrt orig_slice
									#($feature->start + $sub_start),
									($feature->end   + $sub_start),
									$feature->strand,
									$feature->scores,
									);

			#We can have problems here if the original score type
			#does not match the collected score type
			#For max magnitude this is not an issue
			#as we take the larget value from the bin
			#But for other methods this may not be true
			#e.g. count
			#Hence, if we want to preserve the 0 window
			#We must account for this in the feature collector
			#e.g. set_collection_defs_by_ResultSet_window_size?
			#Just omit 0 window for reads

		  }
		}
		
		print "Window size 0 (natural resolution) has ".scalar(@{$features})." feature bins for:\t".$slice->name."\n";	
	  }

	  #Now store bins
	  #	  my ($bin_end, $bin_scores);
	  my $num_bins;

	  foreach my $wsize(sort keys %{$bins}){
		$num_bins = scalar(@{$bins->{$wsize}});
		#warn "$num_bins bin scores for $wsize:\t".join(',', @{$bins->{$wsize}});

		#Should we handle bin trimming here for overhanging slices
		#Then counts wil be correct and wont have to do in caller


		$counts{$wsize}+= $num_bins;	
		
	   

		#We don't need this loop for collections as we can simply push all the scores at once
		#Just use the slice start and end
		if($bin_model eq 'SIMPLE'){

		  $self->write_collection($wsize,
								  $orig_slice,
								  #$sub_start,
								  $sub_end,
								  $orig_slice->strand,#This is most likely 1!
								  #Override this woth 0 in descendant Collector if required.
								  $bins->{$wsize},
								 );

		}
		else{
		  throw('Bio::EnsEMBL::Funcgen::Collector does not yet support non-SIMPLE bin models');
		  #i.e. More than one score
		}
	  

		
#		#Reset start and end for new wsize
#		$bin_start = $slice->start;
#		$bin_end   = $slice->start;
#
#
#
#		#We don't need this loop for collections as we can simply push all the scores at once
#		
#
#		foreach my $bin_index(0..$#{$bins->{$wsize}}){
#
#
#
#		  #default method to handle simple fixed width bin?
#		  #bin_end need to be defined dependant on the bin type
#		  #($bin_start) = $self->process_default_bin($bins->{$wsize}->[$bin_index], $wsize);#?
#
#		  
#
#		  #either define default bin method in descendant
#		  #Or can we set a process_bin_method var?
#		  #No just pass all this info to write collection and handle it there?
#		  
#		  #Can we have just predefined rotueines handling different bin types? 
#		  #Simple
#		  #Simple compressed
#		  #Clipped
#		  #This will prevent hanving to make attrs/method for storing persistent start/end/score info
#		  
#
#
#		  #Need validate bin_type method
#		  #Could convert these to numbers for speed as with binning methods
#
#		  if($bin_model eq 'SIMPLE'){
#
#			$bin_scores = $bins->{$wsize}->[$bin_index];
#
#			warn "bin scores is $bin_scores";
#
#			
#			#next if ! $bin_score;#No we're no inc'ing the start ends for bins with no scores
#
#			$bin_end += $wsize;
#			
#			#if($bin_score){#Removed this as we always want to write the score even if it is 0
#			  
#			  #This is a little backwards as we are generating the object to store it
#			  #If we are aiming for speed the maybe we could also commodotise the store method
#			  #store by args arrays? store_fast?
#			  #Speed not essential for storing!
#			  
#			  #Note: list ref passed
#			  			
#			  #Don't need to pass all this info for fixed width blob collections
#			  #Need to write some default handlers depedant on the collection type
#			  #Simple(original)
#			  #Simple compressed
#			  #Multi compressed
#			  #Clipped uncompressed?
#
#		
#			  $self->write_collection($wsize,
#									  $orig_slice,
#									  ($bin_start + $slice_adj), 
#									  ($bin_end   + $slice_adj),
#									  $orig_slice->strand,#This is most likely 0
#									  $bin_scores,
#									 );
#
#			  #Only count if we have a stored(projected?) feature
#			$count++;#Change this to attr/method?
#			#}
#			
#			$bin_start += $wsize;
#		  }
#		  else{
#			throw('Bio::EnsEMBL::Funcgen::Collector does not yet support non-SIMPLE bin models');
		#	  }
		#	}
		
		#warn "Window size $wsize has ".scalar(@{$bins->{$wsize}})." bins";
		#$counts{$wsize}+= $count;	
	  }
	}

	$store_natural = 0;	#Turn off storing of natural resolution for next chunk length sets
  }
  
  #Now need to write last collections for each wsize
  
  foreach my $wsize(@{$self->window_sizes}){

	next if $wsize == 0 && ! $store_natural;
	next if $wsize != 0 && $only_natural;

	print "Writing final $wsize window_size collection, this may result in slightly different bin numbers from counts due to removing overhanging bins past end of slice\n";

	$self->write_collection($wsize, $orig_slice);#store last collection
  }


  #Print some counts here 
  foreach my $wsize(sort (keys %counts)){
	print "Generated ".$counts{$wsize}." bins for window size $wsize for ".$orig_slice->name."\n";
	#Some may have failed to store if we are projecting to a new assembly
	#Need collection count here too, but would need methods for this?
  }

  #Return this counts hash so we can print/log from the caller, hence we don't print in here?
  
  return;
}



=head2 _bin_features_by_window_sizes

  Args[0]    : Bio::EnsEMBL::Slice
  Args[1]    : ARRAYREF of window sizes
  Args[2]    : int - bin method, currently defined by validate_bin_methods
  Args[3]    : ARRAYREF of Bio::EnsEMBL::Features
  Example    : $bins = $self->_bin_features_by_window_sizes(
													 -slice         => $slice,
													 -window_sizes  => $chunk_windows{$chunk_length},
													 -bin_method    => $bin_method,
													 -features      => $features,
													);
  Description: Bins feature scores for a given list of window sizes and predefined method number
  Returntype : HASHREF of scores per bin per window size
  Exceptions : Throws if bin method not supported
  Caller     : store_window_bins_by_Slice
  Status     : At Risk

=cut


#To do
# 1 Remove Bio::EnsEMBL::Feature dependancy? Or just create Features for non adaptor Collectors.
#   Is there a way we can skip the object generation in the adaptor completely and just 
#   pass the values we need?
# 2 Separate methods, so we can define custom methods in descendants?
# 3 Expand %bins model to optionally be one of
#   the following dependant on binning method
#   Simple:  fixed width containing arrays of scores for each window
#   Multi:   fixed width containing multiple arrays of scores for each window
#   Non-simple?: Separate aggregated features, either fixed width or not, not BLOB!
#   Clipped: default fixed width with option to clip start and end.  Needs start/end attrs
#         Can't store this in a blob due to non-standard start ends?
#         Most likely want more than one score here? Count/Density SNPs?
#         Removes data skew from standard window bins, would need to store each bin and post
#         process. Or do in line to avoid 2nd post-processing loop,requires awareness of when 
#         we have moved to a new bin between features.  This holds for overlapping and 
#         non-overlapping features. Once we have observed a gap we need to clip the end of the
#         last bin and clip the start of the new bin. This requires knowing the greatest end 
#         values from the last bin's feature. what if two overlapping features had the same 
#         start and different end, would we see the longest last? Check default slice_fetch sort

sub _bin_features_by_window_sizes{
  my $this = shift;
  my ( $slice, $window_sizes, $method, $features ) =
    rearrange( [ 'SLICE', 'WINDOW_SIZES', 'BIN_METHOD', 'FEATURES' ], @_ );

  
  #Do this conditional on the Collection type
  #i.e. is collection seq_region blob then no else yes
  #if ( !defined($features) || !@{$features} ) { return {} }

  #warn 'Processing '.scalar(@$features).' features for window sizes '.join(', ',@$window_sizes).' for slice '.$slice->name."\n";	 
  
  #Set up some hashes to store data by window_size
  my (%bins, %nbins, %bin_counts);
  my $slice_start = $slice->start();

   #Default handlers for
  #my($first_bin);
  #if ( $method == 0 ||    # 'count' or 'density'
  #     $method == 3 ||    # 'fractional_count' or 'weight'
  #     $method == 4       # 'coverage'
  #  ){
  #  # For binning methods where each bin contain numerical values.
  #	$first_bin = 0;
  #  } 
  #  else {
  #	# For binning methods where each bin does not contain numerical
  #   # values.
  #
  #	#Remove this
  #	$first_bin = undef;
  #  }


  #Set up some bin data for the windows
  my $slice_length = $slice->length;

  foreach my $wsize (@$window_sizes) {
	#TO DO: Need to modify this block if default 0's are undesirable for collection type
	#i.e. should it be undef instead? May have prolbems representing undef in blob

	$nbins{$wsize}         = int($slice_length / $wsize); #int rounds down
	#nbins is actually the index of the bin not the 'number'
	#Unless slice_Length is a multiple!
	$nbins{$wsize}-- if(! ($slice_length % $wsize));

	#Create default bins with 0
	@{$bins{$wsize}} = ();
	map {$bins{$wsize}->[$_] = 0} (0 .. $nbins{$wsize});		
	
	#Set bin counts to 0 for each bin
	@{$bin_counts{$wsize}}    = ();

	#This is adding an undef to the start of the array!?
    map { $bin_counts{$wsize}->[($_)] = 0 } @{$bins{$wsize}};

	foreach my $bin(@{$bins{$wsize}}){
	  $bin_counts{$wsize}->[$bin] = 0;
	}	
  }

  #warn "bin_counts are :\n".Data::Dumper::Dumper(\%bin_counts);
  #This fails for slices which are smaller than the chunk length; 
  my $feature_index = 0;
  my ($bin_index, @bin_masks);
 
  foreach my $feature ( @{$features} ) {
	#Set up the bins for each window size

	#Omit test for Bio::EnsEMBL::Feature here for speed
	#Only needs start/end methods

	foreach my $wsize (@$window_sizes) {
	
	  #We have already highjacked the object creation by here
	  #This is done in core BaseFeatureAdaptor
	  #We probably don't want to do this for ResultFeatures as we don't use the
	  #standard feature implementation
	  #we already use an array and we don't store the slice
	  #as this is already known by the caller
	  #and we always build on top level so we don't need to remap

	  #We do however need the slice to store, as we only store local starts when generating
	  #We need a store by Slice method?
	  #This will remove the need to inherit from Feature.
	  #These will need to be regenerated everytime we import a new build
	  #As we do with the probe_features themselves
	  #This also mean the result_feature status has to be associated with a coord_system_id
	  
	  #Which bins do the start and end lie in for this feature?
	  #Already dealing with local starts, so no slice subtraction
	  #Could wrap these start/end methods via the descendant Collector
	  #to remove the Feature dependancy? Or just create Features when parsing in the caller
	  my $start_bin =  int(($feature->start ) / $wsize);
	  my $end_bin   =  int(($feature->end) / $wsize );
   	  $end_bin = $nbins{$wsize} if $end_bin > $nbins{$wsize};



	  #Slightly obfuscated code to match method number(faster)
	  #by avoiding string comparisons.
	  #Could call methods directly using coderef set in validate_bin_method
	  #Accessor may slow things down, but should be uniform for all methods
	  #rather than being dependant on position in if/else block below

	  #reserve 0 for descendant defined method?
	  #There fore always fastest in this block, or use coderefs?
	  if ( $method == 0 ) {
		# ----------------------------------------------------------------
		# For 'count' and 'density'.
		
		for ( $bin_index = $start_bin ;
			  $bin_index <= $end_bin ;
			  ++$bin_index ) {

		  $bins{$wsize}->[$bin_index]++;

		  #warn "setting $wsize bin $bin_index to ". $bins{$wsize}->[$bin_index];

		}
	  }

=pod
	
	  } elsif ( $method == 1 ) {
		# ----------------------------------------------------------------
		# For 'indices' and 'index'

	  
		#How is this useful?
		#Is this not just count per bin?
		#No this is a list of the feature indices
		#So forms a distribution?

		throw('Not implemented for method for index'); 

		for ( my $bin_index = $start_bin ;
			  $bin_index <= $end_bin ;
			  ++$bin_index ) {
		  push( @{ $bins[$bin_index] }, $feature_index );
		}

		++$feature_index;

	  } elsif ( $method == 2 ) {
		# ----------------------------------------------------------------
		# For 'features' and 'feature'.

		throw('Not implemented for method for feature/features'); 
		
		for ( my $bin_index = $start_bin ;
			  $bin_index <= $end_bin ;
			  ++$bin_index ) {
		  push( @{ $bins[$bin_index] }, $feature );
		}
		
	  } elsif ( $method == 3 ) {
		# ----------------------------------------------------------------
		# For 'fractional_count' and 'weight'.
		

		throw('Not implemented for method for fractional_count/weight'); 
		
		if ( $start_bin == $end_bin ) {
		  ++$bins[$start_bin];
		} else {
		  
		  my $feature_length =
			$feature->[FEATURE_END] - $feature->[FEATURE_START] + 1;
		  
		  # The first bin...
		  $bins[$start_bin] +=
			( ( $start_bin + 1 )*$bin_length -
			  ( $feature->[FEATURE_START] - $slice_start ) )/
				$feature_length;

		  # The intermediate bins (if there are any)...
		  for ( my $bin_index = $start_bin + 1 ;
				$bin_index <= $end_bin - 1 ;
				++$bin_index ) {
			$bins[$bin_index] += $bin_length/$feature_length;
		  }
		  
		  # The last bin...
		  $bins[$end_bin] +=
			( ( $feature->[FEATURE_END] - $slice_start ) -
			  $end_bin*$bin_length +
			  1 )/$feature_length;
		  
		}						## end else [ if ( $start_bin == $end_bin)

	  } 
	  elsif ( $method == 4 ) {
		# ----------------------------------------------------------------
		# For 'coverage'.
		
		#What exactly is this doing?
		#This is coverage of bin
		#Rather than coverage of feature as in fractional_count
	  
		
	#	my $feature_start = $feature->[FEATURE_START] - $slice_start;
	#	my $feature_end   = $feature->[FEATURE_END] - $slice_start;
	#	
  	#	if ( !defined( $bin_masks[$start_bin] )
	# 		 || ( defined( $bin_masks[$start_bin] )
	#	 		  && $bin_masks[$start_bin] != 1 ) ) {
 	#	  # Mask the $start_bin from the start of the feature to the end
 	#	  # of the bin, or to the end of the feature (whichever occurs
	#	  # first).
 	#	  my $bin_start = int( $start_bin*$bin_length );
 	#	  my $bin_end = int( ( $start_bin + 1 )*$bin_length - 1 ); 
	#	  for ( my $pos = $feature_start;
 	#			$pos <= $bin_end && $pos <= $feature_end ;
 	#			++$pos ) {
 	#		$bin_masks[$start_bin][ $pos - $bin_start ] = 1;
 	#	  }
 	#	}
 	#	
 	#	for ( my $bin_index = $start_bin + 1 ;
 	#		  $bin_index <= $end_bin - 1 ;
 	#		  ++$bin_index ) {
 	#	  # Mark the middle bins between $start_bin and $end_bin as fully
 	#	  # masked out.
 	#	  $bin_masks[$bin_index] = 1;
 	#	}
 	#	
 	#	if ( $end_bin != $start_bin ) {
 	#	
 	#	  if ( !defined( $bin_masks[$end_bin] )
 	#		   || ( defined( $bin_masks[$end_bin] )
 	#				&& $bin_masks[$end_bin] != 1 ) ) {
 	#		# Mask the $end_bin from the start of the bin to the end of
 	#		# the feature, or to the end of the bin (whichever occurs
 	#		# first).
 	#		my $bin_start = int( $end_bin*$bin_length );
 	#		my $bin_end = int( ( $end_bin + 1 )*$bin_length - 1 );
 	#		for ( my $pos = $bin_start ;
 	#			  $pos <= $feature_end && $pos <= $bin_end ;
 	#			  ++$pos ) {
 	#		  $bin_masks[$end_bin][ $pos - $bin_start ] = 1;
 	#		}
 	#	  }
  	#	}
  	 # }							## end elsif ( $method == 4 )

=cut

	  
	  elsif ( $method == 5 ) {
		#$self->$method($bin_index, $start_bin, $end_bin, $wsize, \%bins, \%bin_counts);


		#average score
		#This is simple an average of all the scores for features which overlap this bin
		#No weighting with respect to the bin or the feature
		
		for ( $bin_index = $start_bin ;
			  $bin_index <= $end_bin ;
			  ++$bin_index ) {

		  #we should really push onto array here so we can have median or mean.
		  $bins{$wsize}->[$bin_index] += $this->get_score_by_Feature($feature);
		  $bin_counts{$wsize}->[$bin_index]++;
		}
	  }  
	  elsif( $method == 6){
		#Max magnitude
		#Take the highest value +ve or -ve score
		for ( $bin_index = $start_bin ;
			  $bin_index <= $end_bin ;
			  ++$bin_index ) {

		  #we really need to capture the lowest -ve and higest +ve scores here and post process
		  #To pick between them
		  
		  my $score = $this->get_score_by_Feature($feature);
		  #Write score method as wrapper to scores?

		  $bins{$wsize}->[$bin_index] ||= [0,0]; #-ve, +ve
		  

		  #warn "Comparing wsize $wsize bin $bin_index score $score to ".  $bins{$wsize}->[$bin_index]->[0].' '.$bins{$wsize}->[$bin_index]->[1]."\n";

		  if($score <  $bins{$wsize}->[$bin_index]->[0]){
			#warn "setting -ve bin to $score\n";
			$bins{$wsize}->[$bin_index]->[0] = $score;
		  }
		  elsif($score > $bins{$wsize}->[$bin_index][1]){
			#warn "setting +ve bin to $score\n";
			$bins{$wsize}->[$bin_index]->[1] = $score;
		  }
		}
	  }
	  else {
		throw("Only accomodates average score method");
	  }
	
	
	}
				
  }	## end foreach my $feature ( @{$features...


  #Now do post processing of bins

=pod

  if ( $method == 4 ) {

    # ------------------------------------------------------------------
    # For the 'coverage' method: Finish up by going through @bin_masks
    # and sum up the arrays.

    for ( my $bin_index = 0 ; $bin_index < $nbins ; ++$bin_index ) {
       if ( defined( $bin_masks[$bin_index] ) ) {
        if ( !ref( $bin_masks[$bin_index] ) ) {
          $bins[$bin_index] = 1;
        } else {
          $bins[$bin_index] =
            scalar( grep ( defined($_), @{ $bin_masks[$bin_index] } ) )/
            $bin_length;
        }
      }
    }
  }

=cut

  if( $method == 5){
	#For average score, need to divide bins by bin_counts

	foreach my $wsize(keys %bins){

	  foreach my $bin_index(0..$#{$bins{$wsize}}){

		if($bin_counts{$wsize}->[$bin_index]){
		  $bins{$wsize}->[$bin_index] /= $bin_counts{$wsize}->[$bin_index];
		}
		#warn "bin_index $wsize:$bin_index has score ".$bins{$wsize}->[$bin_index];
	  }
	}
  }
  elsif( $method == 6){
	#Max magnitude
	#Take the highest value +ve or -ve score

	foreach my $wsize(keys %bins){

	  foreach my $bin_index(0..$#{$bins{$wsize}}){

		#So we have the potential that we have no listref in a given bin

		#default value if we haven't seen anything is 0
		#we actually want an array of -ve +ve values

		#warn "Are we storing 0 values for absent data?";
		#Not for max_magnitude, but maybe for others?

		if($bins{$wsize}->[$bin_index]){
		  #warn $wsize.':'.$bin_index.':'.$bins{$wsize}->[$bin_index]->[0].'-'.$bins{$wsize}->[$bin_index]->[1];
		  my $tmp_minus = $bins{$wsize}->[$bin_index]->[0] * -1;
		  
		  if($tmp_minus > $bins{$wsize}->[$bin_index]->[1]){
			$bins{$wsize}->[$bin_index] = $bins{$wsize}->[$bin_index]->[0];
		  }
		  else{
			$bins{$wsize}->[$bin_index] = $bins{$wsize}->[$bin_index]->[1];
		  }

		  #warn "bin $bin_index now ".	$bins{$wsize}->[$bin_index];
		}
	  }
	}
  }
  elsif($method != 0){#Do no post processing for count(0)
	throw('Collector currently only accomodates average_score, count and max magnitude methods');
  }


  #Could return bin_counts too summary reporting in zmenu
  #Could also do counting of specific type

  #warn "returning bins ".Data::Dumper::Dumper(\%bins);

  return \%bins;
} ## end sub _bin_features


=pod

#These could potentially be used as code refs to avoid having the if else block
#This way we can also define new methods in the descendant Collector?
#Would have to have pass args and refs to bin hashes
#This would slow things down over direct access here
#But speed is no longer that critical as we do not use the Collector for display
#purposes, only to build the Collections which are then used for display directly.

sub calculate_average_score{
  my $self = shift;

  if ( $method == 5 ) {
		#average score
		#This is simple an average of all the scores for features which overlap this bin
		#No weighting with respect to the bin or the feature
		
		for ( my $bin_index = $start_bin ;
			  $bin_index <= $end_bin ;
			  ++$bin_index ) {

		  #we should really push onto array here so we can have median or mean.
		  $bins{$wsize}->[$bin_index] += $feature->score;
		  $bin_counts{$wsize}->[$bin_index]++;
		}
	  }  



}
	

sub post_process_average_score{

}


sub calculate_max_magnitude{
  my $self = shift;

  #Max magnitude
  #Take the highest value +ve or -ve score
  for ( my $bin_index = $start_bin ;
		$bin_index <= $end_bin ;
		++$bin_index ) {
	
	#we really need to capture the lowest -ve and higest +ve scores here and post process
	#To pick between them
	
	my $score = $feature->score;
	$bins{$wsize}->[$bin_index] ||= [0,0]; #-ve, +ve
	
	if($score <  $bins{$wsize}->[$bin_index]->[0]){
	  $bins{$wsize}->[$bin_index]->[0] = $score;
	}
	elsif($score > $bins{$wsize}->[$bin_index][1]){
	  $bins{$wsize}->[$bin_index]->[1] = $score;
	}
  }
}


sub post_process_max_magnitude{

}

=cut

#separated to allow addition of non-standard methods
#Could potentially add these in new
#and put this back in _bin_features


sub validate_bin_method{
  my ($self, $method) = @_;


  #change this to set the coderefs
  #Just set anonymous sub to immediately return for non post processed methods
  #No need for coderef, just set the method name?

  #if(! $self->can('calculate_'.$method)){
  #throw("$method method does not have a valid calculate_${method} method");
  #}

  #if($self->can('post_process_'.$method)){
  ##set post process flag?
  #or simply do this can in line in the _bin_features sub?
  #}
  



  #Add average_score to avoid changing Collection.pm
  my $class = ref($self);
  ${$class::VALID_BINNING_METHODS}{'average_score'} = 5;
  ${$class::VALID_BINNING_METHODS}{'max_magnitude'} = 6;
  ${$class::VALID_BINNING_METHODS}{'count'} = 0;
  


    #foreach my $method_name(keys %{$class::VALID_BINNING_METHODS}){
#	warn "valid method is $method name";
#  }


  if ( ! exists( ${$class::VALID_BINNING_METHODS}{$method} ) ) {
    throw(
		  sprintf(
				  "Invalid binning method '%s', valid methods are:\n\t%s\n",
				  $method,
				  join( "\n\t", sort( keys(%{$class::VALID_BINNING_METHODS}) ) ) ) );
  }
  else{
	#warn "found valid method $method with index ".${$class::VALID_BINNING_METHODS}{$method};
  }
  
  return ${$class::VALID_BINNING_METHODS}{$method};
}



1;
