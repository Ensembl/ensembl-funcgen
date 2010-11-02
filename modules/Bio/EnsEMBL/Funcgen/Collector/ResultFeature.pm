# $Id: ResultFeature.pm,v 1.8.8.1 2010-11-02 09:52:21 nj1 Exp $



=head1 LICENSE

  Copyright (c) 1999-2009 The European Bioinformatics Institute and
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


package Bio::EnsEMBL::Funcgen::Collector::ResultFeature;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument  ('rearrange');
use Bio::EnsEMBL::Utils::Exception ('throw');
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw (open_file);
use Bio::EnsEMBL::Funcgen::Collection::ResultFeature;

use base('Bio::EnsEMBL::Utils::Collector');#ISA


### Global config variables
# See associated Collector methods for more info

$Bio::EnsEMBL::Utils::Collector::bin_model  = 'SIMPLE';

#Default is for read coverage, array intensity config defined in set_collection_defs_by_ResultSet
#Kept here as example of basic package config
$Bio::EnsEMBL::Utils::Collector::window_sizes = [30, 65, 130, 260, 450, 648, 950, 1296];


#May need to drop 30 here due to RPKM float value doubling packed size

# Tuning 900(Most used display size) - 172(non-drawable area) = 772 pixels/bins
# Optimal window sizes over zoom levels
# 1kb 5kb 10kb 50kb 100kb 200kb 500kb 1mb
# 2   7   13   65   130   260   648   1296
#Theoretical max for default 16MB max_allowed_packet_size is 30(~23kb) - sufficient min for read coverage.
#0 will be used for array intensies until 65 is more appropriate(see adaptor)
#450 and 950 added to handle middle ground, but are these really needed?
#New windows will use more memory due to smaller window sizes
#But we have dropped a window size from the mid-upper range
#Largest size is essential so we aren't retrieving/drawing more features than pixels
#Smaller sizes are desirable as more resolution required when zoomed in.

$Bio::EnsEMBL::Utils::Collector::bin_method   = 'RPKM';#'count';#only used by collector
#small integers(little endian)
#RPKM is now float(native)
$Bio::EnsEMBL::Utils::Collector::packed_size   = 4;#2;#per score
$Bio::EnsEMBL::Utils::Collector::pack_template = 'f';#'v';#per score


### Mandatory methods required by the base Collector




# You must have a methods like this to perform the store

sub store_window_bins_by_Slice_ResultSet {
  my ($self, $slice, $rset, %config) = @_;

  $self->source_set_type('result');#required by get_Feature_by_Slice   
  $self->set_collection_defs_by_ResultSet($rset);  
  $self->set_config(%config);

  $self->store_window_bins_by_Slice($slice);
  
  return;
}



sub store_window_bins_by_Slice_Parser{
  my ($self, $slice, $imp, %config) = @_;


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
  
  my ($skip_zero_window) = rearrange( [ 'SKIP_ZERO_WINDOW'], %config );



  #There is no way of setting IMPORTED status for slice based jobs
  #We need the accumulator to set the IMPORTED/RESULT_FEATURE_SET status
  #Would need to turn parse_and_import into RunnableDB?
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
 
  $self->set_config(%config, (-method_config => {
												 #RPKM method config
												 -dnadb          => $imp->db->dnadb,
												 -total_features => $imp->total_features,
												 -gender         => $imp->cell_type->gender,
												}
							 ));

  $self->store_window_bins_by_Slice($slice);
  
  return;
}


=head2 rollback_Features_by_Slice

  Args[0]    : Bio::EnsEMBL::Slice
  Example    : $collector->get_Feature_by_Slice($slice);
  Description: Wrapper method to fetch input features for building the Collections
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::Collection::ResultFeatures
  Exceptions : None
  Caller     : Bio::EnsEMBL::Utils::Collector::store_window_bins_by_Slice
  Status     : At Risk

=cut

#optional method to perform feature rollback

sub rollback_Features_by_Slice{
  my ($self, $slice) = @_;

  #Point to Helper here
  #This is already done in the InputSet importer for
  #bed/sam seq imports
  #but not for array based imports
  #Need to take account of wsizes


}



=head2 get_Features_by_Slice

  Args[0]    : Bio::EnsEMBL::Slice
  Example    : $collector->get_Feature_by_Slice($slice);
  Description: Wrapper method to fetch input features for building the Collections
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::Collection::ResultFeatures
  Exceptions : None
  Caller     : Bio::EnsEMBL::Utils::Collector::store_window_bins_by_Slice
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
  Caller     : Bio::EnsEMBL::Utils::Collector::store_window_bins_by_Slice
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
	
	
	if($new_cps >= $self->max_data_type_size){
	  warn("Have found $wsize collection larger that max_data_type_size(16MB). Need to implement cross collection querying");
	  #Via a union of two substr queries
	  #This is no set to 64MB, but this does not directly translate to the maximum size allowed in a single insert
	}
	  


	if($wsize == 0){
	  #We need to 0 wsize collections to be projected as single features
	  #Other wise the assembly projection will not work
	  #Then use the 0 wsize ResultFeatures to generate

	  #if($new_cps >= $self->max_data_type_size){
	  #	warn('Have found collection larger that max_data_type_size(16MB). Need to implement cross collection querying');
	  #	#Via a union of two substr queries
	  #  }
	  

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

=head2 get_score_by_Feature


=cut 

#Not required for count/density based collections

sub get_score_by_Feature{
  my ($self, $feature) = @_;

  #For speed, assume we have a valid Bio::EnsEMBL::Funcgen::Collection::ResultFeature
  #This will always be a 0 wsize collection, hence we always want the first score

  return $feature->scores->[0];

}



=head2 reinitialise_input 


=cut 

#optional method - resets file handle for iterative slice parsing for different chunk length sets

sub reinitialise_input{
  my $self = shift;

  #No need to do this for DB queries
  return if $self->source_set_type eq 'result';

  #This will make all the parser counts go screwy
  #as we are potentially counting the same features again


  #Move this to parser method
  $self->parser->file_handle->close;
  $self->parser->file_handle( open_file($self->parser->input_file, $self->parser->input_file_operator) );

  #Reset caches 
  $self->parser->{'last_slice'} = undef;
  $self->parser->{'overhang_features'} = [];
  $self->parser->last_line('');
  return;
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
  Caller     : write_collection and Bio::EnsEMBL::Utils::Collector::store_window_bins_by_Slice
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
  Example    : $self->set_collection_defs_by_ResultSet($rset);
  Description: Sets the packed_size and pack_template and bin_method dependant on the 
               ResultSet type(reads/array)
  Returntype : None
  Exceptions : throws if supporting InputSet is not of type result (i.e. short reads import)
               throws if supporting InputSet format is not SEQUENCING (i.e. short reads import)
               throws if ResultSet is not and input_set or experimental_chip based ResultSet (i.e. channel etc)
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
	
	if($rset->table_name eq 'experimental_chip'){ #Array Intensities i.e. single float
	  #perl only offers native endian order for floats (http://www.perlmonks.org/?node_id=629530)
	  #$Bio::EnsEMBL::Utils::Collector::packed_size       = 4;#per score
	  #$Bio::EnsEMBL::Utils::Collector::pack_template     = 'f';#per score
	  #Use defaults for these now

	  $Bio::EnsEMBL::Utils::Collector::bin_method        = 'max_magnitude';#only used by collector
	  $Bio::EnsEMBL::Utils::Collector::window_sizes->[0] = 0;#Can have natural resolution for low density array data
	}
	elsif($rset->table_name eq 'input_set'){
	  #Need to reset this as we may be doing serial queries.
	  $Bio::EnsEMBL::Utils::Collector::window_sizes->[0] = 30;

	  #Currently only expecting int from InputSet
	  my @isets = @{$rset->get_InputSets};
	  my @tmp_isets = grep(!/result/, (map $_->feature_class, @isets));
	  
	  if(@tmp_isets){
		throw("Bio::EnsEMBL::Funcgen::Collector::ResultFeature only supports result type InputSets, not @tmp_isets types");
	  }

	  #We still have no way of encoding pack_type for result_feature InputSets

	  @tmp_isets = grep(!/SEQUENCING/, (map $_->format, @isets));
	  
	  if(@tmp_isets){
		throw("Bio::EnsEMBL::Funcgen::Collector::ResultFeature only supports SEQUENCING format InputSets, not @tmp_isets formats");
	  }
	}
	else{
	  throw('Bio::EnsEMBL::Funcgen::Collector:ResultFeature does not support ResultSets of type'.$rset->table_name);
	}
  }

  #Do we need to validate the smallest non-0 window size
  #against the max pack size?
  #This should be done in the Collector

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
