# $Id: ResultFeature.pm,v 1.5 2010-03-18 10:21:28 nj1 Exp $

package Bio::EnsEMBL::Funcgen::Collector::ResultFeature;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument  ('rearrange');
use Bio::EnsEMBL::Utils::Exception ('throw');
use Bio::EnsEMBL::Funcgen::Collection::ResultFeature;

use base( 'Bio::EnsEMBL::Funcgen::Collector');#ISA

#Removed these as we now set_collection_def_by_ResultSet
#$Bio::EnsEMBL::Funcgen::Collector::intensity_bin_method = 'max_magnitude';
#$Bio::EnsEMBL::Funcgen::Collector::read_bin_method = 'count';

#Read Coverage i.e. small integers(little endian)
#$Bio::EnsEMBL::Funcgen::Collector::read_packed_size   = 2;#per score
#$Bio::EnsEMBL::Funcgen::Collector::read_pack_template = 'v';#per score

#Array Intensities i.e. single float(native)
#$Bio::EnsEMBL::Funcgen::Collector::intensity_packed_size = 4;#per score
#$Bio::EnsEMBL::Funcgen::Collector::intensity_template    = 'f';#per score
#This is only true for int values i.e. read coverage
#For array value we need a single float
#perl only offers native endian order for floats
#Good summary here about this:
#http://www.perlmonks.org/?node_id=629530



$Bio::EnsEMBL::Funcgen::Collector::bin_model  = 'SIMPLE';
#Only the packed size is required in the collector
#Packing done in storing object(e.g.adaptor)
#Define here if this is the storing object
#Otherwise in the adaptor as package cars or methods to override those in base Collector
#$super_class::pack_template = 'v';#per score #Move to Collection::ResultFeature?
#$super_class::packed_size   = '2';#per score
#$super_class::window_sizes  = [0, 150, 250, 350, 450, 550, 650];#Can remove this from ResultFeatureAdaptor?
#Window sizes are actually used in the adaptor, so another case
#where we have an adaptor override. Only specify this if you are not using an adaptor




#Can we write default SimpleCollection objects etc
#So we don't need to write a collection for uncompressed/collected data
#This will need generic methods in base Collector?
#Would also need to put default config in adaptor

#Using standard object generation methods for ResultFeature
#We never want the heavy object from the result_feature table method
#In fact we never want the heavy object unless we create it from scratch to store it
#Therefore we don't need to use create_feature at all unless we ever want to use
#The collection to generate bins on the fly(rather than for storing)
#We would only want this if we don't want to implement the object, just the array
#Now we're back to considering the base collection array
#'DBID', 'START', 'END', 'STRAND', 'SLICE'
#we don't need dbID or slice for ResultFeature
#But we do need a score
#Stick with object implementation for now i.e. dont use create_feature methods


#TO DO
#Had to change ResultFeature to inherit from Feature, so can use Feature methods
#Is this going to cause problem mixing and array based object with hash based object?
#Not is we implement the mandatory methods here.
#seq_region start etc?

#Issues
#As we define the windows by the start of the features
#we get some skew in the data to the right as we generally have 
#overhanging 3' features and no overhanging 5' features
#Also when dealing with window sets we are not resetting the start 
#of the bins when we reach a gap, hence we will then get 5' skew
#but only for gaps within a chunk, not at the start.
#Solution, set start and end values for each bin dynamically
#to the start or end of a feature if there are not boundary spanning features
#This should remove skew completely


### These first methods are used in the Collector and are required unless otherwise stated

# You must have a method like this to perform the store

sub store_window_bins_by_Slice_ResultSet {
  my ($self, $slice, $rset, %config) = @_;

  #Need to 'set_config' in caller
  #Or should we do this once in the caller
  #Then we can do tests here before calling super method

  $self->source_set_type('result');#required by get_Feature_by_Slice   
  $self->set_collection_defs_by_ResultSet($rset);
 
  #When called this does not pass $self, so we would have to pass the result set explicitly
  #Let's not pass a code ref, let's just write a wrapper instead
  #my $fetch_method_ref = $rset->can('get_ResultFeatures_by_Slice');
  #We could pass a ref to the wrapper method? But this is unlikely to enable direct usage 
  #of 3rd party parser method
  
  $self->store_window_bins_by_Slice($slice, %config);#, (
  
  return;
}

#This will be called by the Importer
#From the Bed parser. Which will contain the methods to access the features
#Hence will need to pass the Bed parser itself to provide access to the parse methods

#This is really more like store_window_bins_by_Slice_Importer
#parser should have InputSet and ResultSet defined

sub store_window_bins_by_Slice_Parser{
  my ($self, $slice, $imp, %config) = @_;

  #Other params?
  #imp, which has had filehandle, ResultSet and InputSet already established
  #Need to bsub these slice jobs, so parse_and_import.pl needs to handle submitting to farm
  #We also need to define window_sizes = 0
  #If we are loading data on an old assembly
  #Then run populate_result_features?
  
  #This depends on the density of the data!
  #For ChIP-Seq alignments we don't want to load all the raw reads
  #Hence we don't want 0 window_size
  #and will need all the input remapped to the required assembly before
  #we import here

  #We need to test for new_assm if we have a reads result set and fail
  #currently don't have access to new assm for validation?
  #Can project to new_assm as this would require two passes storing initially on the 0 window level
  #which we don't want to do
  #throw if new_assm as we need to remap before running this
  #This is done in Collector so why are we parsing it here?
  #skip_zero window will always be set, so this would fail!!!

  #test parse config here, will this strip it out of the hash before
  #we pass it to the super method?
  #Also need to test skip_zero_window or window_sizes dependant on result set type(sequencing, array)
  #This needs to be defined when craeting the ResultSet in the Importer
  #NEED TO CHANGE ALL ResultSet generation to add type!!!!!
  
  my ($wsizes, $new_assm, $skip_zero_window) =
    rearrange( [ 'WINDOW_SIZES', 'NEW_ASSEMBLY', 'SKIP_ZERO_WINDOW'], %config );


  #Can/should we think about pipelining this?
  #For a ChIP Seq set this would be two or three jobs with varying interdependencies
  #1 Align to genome (bwa/maq/bowtie)
  #2 Load ResultFeatures (depends on 1)
  #3 Call peaks (depends on 1)
  
  #There is no way of setting IMPORTED status for slice based jobs
  #As we will need to wait for all to finish
  
  #We need the accumulator to set the IMPORTED status



  #Would need to turn parse_and_import into RunnableDB
  #To run these as slice jobs we would need to code to configure some slice based input ids
  

 
  $self->source_set_type('input');#required by get_Feature_by_Slice
  $self->set_collection_defs_by_ResultSet($imp->result_set);
  $self->parser($imp);
 
  #For safety, set skip_zero window if we are using SEQUENCING data
    
  if(! $skip_zero_window){
	#Assume we only have one set here(enforced in define_and_validate_sets)
	my ($iset) = @{$imp->result_set->get_InputSets};

	if($iset->format eq 'SEQUENCING'){
	  $config{'-skip_zero_window'} = 1;
	}
  }
 

  #When called this does not pass $self, so we would have to pass the result set explicitly
  #Let's not pass a code ref, let's just write a wrapper instead
  #my $fetch_method_ref = $rset->can('get_ResultFeatures_by_Slice');
  #We could pass a ref to the wrapper method? But this is unlikely to enable direct usage 
  #of 3rd party parser method
  
  $self->store_window_bins_by_Slice($slice, %config);
  
  return;
}


=head2 rollback_Features_by_Slice

  Args[0]    : Bio::EnsEMBL::Slice
  Example    : $collector->get_Feature_by_Slice($slice);
  Description: Wrapper method to fetch input features for building the Collections
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::Collection::ResultFeatures
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::Collector::store_window_bins_by_Slice
  Status     : At Risk

=cut

#optional method to perform feature rollback

sub rollback_Features_by_Slice{
  my ($self, $slice) = @_;

  #Point to Helper here
  #This is already done in the InputSet importer for
  #seq imports
  #but not for array based imports
  #Need to take account of wsizes


}



=head2 get_Features_by_Slice

  Args[0]    : Bio::EnsEMBL::Slice
  Example    : $collector->get_Feature_by_Slice($slice);
  Description: Wrapper method to fetch input features for building the Collections
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::Collection::ResultFeatures
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::Collector::store_window_bins_by_Slice
  Status     : At Risk

=cut

#Can remove this method if we pass a method ref
#Default for normal FeatureAdaptors would fetch_all_by_Slice
#However cannot do this as we don't get self passed via a code ref
#we could however manually pass this when we call the method?
#self will either be a collector or a parser
#Just maintain this wrapper for now for simplicity

sub get_Features_by_Slice{
  my ($self, $slice) = @_;

  #Can we pass this as method ref to Collector, this would prvent the need to write
  #this wrapper for standard fetch_all_by_Slice based access
  my $features;
  my $source_set_type = $self->source_set_type;

  if ($source_set_type eq 'result'){
	#Add more args here? status?
	#Add default window size 0, if ResultSet is alread a RESULT_FEATURE_SET
	#This may occur if you are projecting features from an old assembly to 
	#0 window size before generating the other windows in a second pass
	$features = $self->result_set->get_ResultFeatures_by_Slice($slice, undef, undef, undef, 0);
  }
  elsif($source_set_type eq 'input'){
	$features = $self->parser->parse_Features_by_Slice($slice);

	#This method assumes a sorted file handle
	#Can we set markers for disk seeking on a sorted handle?
	#Either we make it handle an slice passed
	#Or we assume the next query slice will be after the last
	#Just restrict to one slice at a time for now
	
	#Also needs to hold cache of long features
	#Kind of reinventing the BaseFeatureAdaptor wheel here?


  }
  else{
	#We should have already validated this by now
	throw('get_Features_by_Slice only support input and result set_types');
  }

  return $features;

}



=head2 write_collection

  Args[0]    : int : window_size
  Args[1]    : Bio::EnsEMBL::Slice
  Args[2]    : optional int : Feature end
  Args[3]    : optional int : Feature strand
  Args[4]    : optional ARRAYREF of scores
  Example    : $self->store_collection($wsize, $slice, $self->collection_start($wsize), $self->collection_end($wsize), $self->collection_strand($wsize));
  Description: Writes Collections by caching scores, defining collection end and strand and 
               storing when end of collection seen or max_packed_size is exceeded (currently does 
               not support fragmented collections). Stores last Collection if only window_size
               and Slice are passed.
  Returntype : None
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::Collector::store_window_bins_by_Slice
  Status     : At Risk

=cut

sub write_collection{
  my ($self, $wsize, $slice, $slice_end, $strand, $scores) = @_;

  if(defined $slice_end){
	#collection_start is defined in Collector for first bin in collection
	$self->collection_end($wsize, $slice_end);
	$self->collection_strand($wsize, $strand);
  }
  

  if(defined $scores){
	my $new_cps = $self->current_packed_size($wsize) + (scalar(@$scores)*$self->packed_size);	
	#This loop is very similar to the _obj_from_sth loop
	#We need to set the slice and modify the start end appropriately
	
	if(($new_cps >= $self->max_data_type_size) ||
	  $wsize == 0){
	  #We need to 0 wsize collections to be projected as single features
	  #Other wise the assembly projection will not work
	  #Then use the 0 wsize ResultFeatures to generate

	  if($new_cps >= $self->max_data_type_size){
		warn('Have found collection larger that max_data_type_size(16MB). Need to implement cross collection querying');
	  }
	  

	  #Need to cache current score for 0 window
	  $self->score_cache($wsize, $scores) if $wsize == 0;
	  $self->store_collection($wsize, $slice, $self->collection_start($wsize), $slice_end, $strand);

	}
	
	$self->score_cache($wsize, $scores) if $wsize != 0; 
  }
  else{#No score info, so we store the remaining last collections for each wsize
	$self->store_collection($wsize, $slice, $self->collection_start($wsize), $self->collection_end($wsize), $self->collection_strand($wsize));
  }

  return;

}


sub get_score_by_Feature{
  my ($self, $feature) = @_;

  #For speed, assume we have a valid Bio::EnsEMBL::Funcgen::Collection::ResultFeature
  #This will always be a 0 wsize collection, hence we always want the first score

  return $feature->scores->[0];

}

##############################################################
### Following methods are not mandatory for a Collector
### but support above methods for this ResultFeature Collector
##############################################################


=head2 store_collection

  Args[0]    : int : window_size
  Args[1]    : Bio::EnsEMBL::Slice
  Args[2]    : int : Feature start
  Args[3]    : int : Feature end
  Args[4]    : int : Feature strand
  Example    : $self->store_collection($wsize, $slice, $self->collection_start($wsize), $self->collection_end($wsize), $self->collection_strand($wsize));
  Description: Collection storage method. Resets seq_regios_Start/end and scores appropriately
               if collection exceeds storage slice. Resets score cache and next collection_start.
  Returntype : None
  Exceptions : None
  Caller     : write_collection and Bio::EnsEMBL::Funcgen::Collector::store_window_bins_by_Slice
  Status     : At Risk

=cut

#Or could print to file

sub store_collection{
  my ($self, $wsize, $full_slice, $slice_start, $slice_end, $strand) =  @_;
	  
  #warn "Storing collection $wsize, $full_slice, $slice_start, $slice_end, $strand";


  my $sr_start = $slice_start;
  my $sr_end   = $slice_end;
  $strand = 0 if ! defined $strand;
  #Overwriting strand value with 0 
  #As we have collected features from both strands?
  $strand = 0 if $wsize != 0;


  #This happens if the last collection was already written in the loop
  #Handle this here so we don't have to set/test collection_start for every record
  return if($sr_start == ($sr_end + 1));


  #Set store slice to start at 1 (PARs, test slices)
  #Can we remove this store_slice now?
  #How are we going to handle storing/fetching on PARs?
  my $store_slice = $full_slice;

  
  if($store_slice->start != 1){
	#Alter sr_start/end
	$sr_start  = $slice_start + $full_slice->start - 1;
	$sr_end = $slice_start + $full_slice->start - 1;
	$store_slice = $store_slice->adaptor->fetch_by_region(undef, $store_slice->seq_region_name);
  }


  ### Splice scores if collection overhang store slice
  #We reassign this below to avoid any weird ref updating behaviour
  my $scores_ref = $self->score_cache($wsize);

  if($sr_end > $full_slice->end){
	#Can happen for 0 wsize if we are running with small test slice
	#But we never want to splice scores for 0 wsize

	if($wsize == 0){  #Throw away any which are not at least partially on this slice
	  return if $sr_start > $full_slice->end;
	}
	else{

	  #Trim the scores and sr_end to the nearest bin
	  my $tmp_end = int($full_slice->end/$wsize) * $wsize;
	  $tmp_end += $wsize if $tmp_end < $full_slice->end;
	  my $overhang_length = ($sr_end - $tmp_end)/$wsize;
	  #This should always be an int as $slice_end is always a valid bin end
	  $sr_end = $tmp_end;
	  
	  #Now remove the surplus scores from the cache
	  splice(@$scores_ref, ($#{$scores_ref} - $overhang_length), $overhang_length);
	}

  }
 
  print 'Storing '.scalar(@$scores_ref)." bins(window_size=$wsize) for:\t".$store_slice->name."\n";

  $self->store([Bio::EnsEMBL::Funcgen::Collection::ResultFeature->new_fast({
																			start  => $sr_start,
																			end    => $sr_end,
																			strand => $strand,
																			scores => $scores_ref,
																			#probe  => undef,
																			result_set_id => $self->result_set->dbID,
																			window_size   => $wsize,
																			slice         => $store_slice,
																		   }
																		  )], $self->result_set, $self->new_assembly);
  
  #Reassign/clean score_cache reference to avoid any reference updating problems
 
  $self->{'score_cache'}{$wsize} = [];

  #Set collection start, this is reset to the actual feature start for wsize == 0
  $self->collection_start($wsize, ($slice_end+1));
 
  return;
}



=head2 set_collection_defs_by_ResultSet

  Args[0]    : Bio::EnsEMBL::Funcgen::ResultSet
  Example    : $self->set_collection_defs($rset);
  Description: Sets the packed_size and pack_template and bin_method dependant on the 
               ResultSet type(reads/array)
  Returntype : None
  Exceptions : None
  Caller     : store_window_bins_by_Slice_ResultSet
               Bio::EnsEMBL::Funcgen:DBSQL::ResultFeatureAdaptor::fetch_all_by_Slice_ResultSet
  Status     : At Risk

=cut

sub set_collection_defs_by_ResultSet{
  my ($self, $rset) = @_;
  
  $self->result_set($rset);

  if(defined $rset){
	#This is most likely already done in the caller
	#$self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::ResultSet', $rset);

	#Keep package vars for clarity?
	
	#This is called by fetch_all_by_Slice_ResultSet
	#So we have to use $rset->table_name instead of source_set_type
	
	if($rset->table_name eq 'experimental_chip'){#i.e. is normal ResultSet with float result
	  $self->{'packed_size'}   = 4;
	  $self->{'pack_template'} = 'f';
	  $self->{'bin_method'}    = 'max_magnitude';#only used by collector
	  $self->{'window_sizes'}  = [0, 150, 300, 450, 600, 750, 900, 1150];
	}
	elsif($rset->table_name eq 'input_set'){
	  #Currently only expecting int from InputSet
	  my @isets = @{$rset->get_InputSets};
	  my @tmp_isets = grep(!/result/, (map $_->feature_class, @isets));
	  
	  if(@tmp_isets){
		throw("Bio::EnsEMBL::Funcgen::Collector::ResultFeature only supports result type InputSets, not @tmp_isets types");
	  }

	  #We still have no way of encoding pack_type for result_feature InputSets
	  #

	  @tmp_isets = grep(!/SEQUENCING/, (map $_->format, @isets));
	  
	  if(@tmp_isets){
		throw("Bio::EnsEMBL::Funcgen::Collector::ResultFeature only supports SEQUENCING format InputSets, not @tmp_isets formats");
	  }

	  $self->{'packed_size'}   = 2;
	  $self->{'pack_template'} = 'v';
	  $self->{'bin_method'}    = 'count';
	  $self->{'window_sizes'}  = [50, 150, 300, 450, 600, 750, 900, 1150];
	}
	else{
	  throw('Bio::EnsEMBL::Funcgen::Collector:ResultFeature does not support ResultSets of type'.$rset->table_name);
	}		
  }

  #Do we need to validate the smallest non-0 window size
  #against the max pack size?
  #This should be don in the Collector

  #warn "Collection defs are:\n".
#	"\tpacked_size:\t". $self->{'packed_size'}."\n".
#	  	"\tpack_template:\t". $self->{'pack_template'}."\n".
#		  "\tbin_method:\t". $self->{'bin_method'}."\n".
#			"\twindow_sizes:\t". join(',', @{$self->{'window_sizes'}})."\n";

  return;

}



=head2 result_set

  Args[0]    : optional Bio::EnsEMBL::Funcgen::ResultSet
  Example    : $self->store_collection($wsize, $slice, $self->collection_start($wsize), $self->collection_end($wsize), $self->collection_strand($wsize));
  Description: result_set attribute method required for fetch wrapper method
  Returntype : Bio::EnsEMBL::Funcgen::ResultSet
  Exceptions : throws if arg is not valid
  Caller     : get_Features_by_Slice
  Status     : At Risk

=cut

sub result_set{
  my ($self, $rset) = @_;

 
  #Can't use is_stored_and_valid here

  if($rset && ! (ref($rset) && $rset->isa('Bio::EnsEMBL::Funcgen::ResultSet'))){
	throw('You must pass a valid Bio::EnsEMBL::Funcgen::ResultSet');
  }
  elsif($rset){

	#Now done in the caller
	#if($rset->has_status('RESULT_FEATURE_SET')){
	#throw('ResultSet('.$rset->name.') already has precomputed ResultFeatures stored, please rollback ResultFeature first');
	  #Retrieving ResultFeatures here would end up retrieving them from the result_feature table
	  #which is where we want to store them
	  #They are only retrived from the probe/result/probe_feature table if they do not have this status.
	  #REMEMBER TO SET THIS IN THE CALLING SCRIPT!?
	#}
  
	$self->{'result_set'} = $rset;
  }
  
  return $self->{'result_set'};
}


=head2 parser

  Args[0]    : optional Bio::EnsEMBL::Funcgen::Parsers::InputSet
  Example    : $self->parser($parser);
  Description: Getter/Setter for parser attribute if this ResultFeature Collector
  Returntype : Bio::EnsEMBL::Funcgen::Parsers::InputSet
  Exceptions : throws if arg is not valid
  Caller     : general
  Status     : At Risk

=cut

sub parser{
  my ($self, $parser) = @_;
 
  #Can't use is_stored_and_valid here

  if($parser && ! (ref($parser) && $parser->isa('Bio::EnsEMBL::Funcgen::Parsers::InputSet'))){
	throw('You must pass a valid Bio::EnsEMBL::Funcgen::Parsers::InputSet');
  }
  elsif($parser){
  
	$self->{'parser'} = $parser;
  }
  
  return $self->{'parser'};
}




sub source_set_type{
  my ($self, $type) = @_;

  $self->{source_set_type} = $type if $type;

  return $self->{source_set_type};
}


1;
