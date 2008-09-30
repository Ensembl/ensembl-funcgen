# $Id: ResultFeature.pm,v 1.6 2008-09-30 10:51:22 nj1 Exp $

package Bio::EnsEMBL::Funcgen::Collection::ResultFeature;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument  ('rearrange');
use Bio::EnsEMBL::Utils::Exception ('throw');
use Bio::EnsEMBL::Funcgen::ResultFeature;

use base( 'Bio::EnsEMBL::Funcgen::DBSQL::ResultFeatureAdaptor',
		  'Bio::EnsEMBL::Collection');

#Had to put adaptor first as it wasn't finding new?
#Wasn't transvering the package tree?



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


sub store_window_bins_by_Slice_ResultSet {
  my $this = shift;
  my ( $slice, $rset, $window_sizes, $logic_name, $rollback) =
    rearrange( [ 'SLICE', 'RESULT_SET', 'WINDOW_SIZES', 'LOGIC_NAME', 'ROLLBACK'], @_ );


  my $slice_adaptor = $this->db->dnadb->get_SliceAdaptor;
  my $method='average_score';#average
  #Hardcoded method to average.  What about those probe which cross boundaries?
  #Will this be returned?

  if(! defined $window_sizes){

	if($this->can('window_sizes')){
	  @{$window_sizes} = $this->window_sizes;
	  warn 'No window sizes provided running with defaults: '.join(', ', @$window_sizes);
	}
	else{
	  throw('No default window sizes available. Must supply at least one window size as an arrayref');
	}
  }
  elsif(scalar(@$window_sizes) == 0){
	throw('Must supply at least one window size as an arrayref, or leave undef for defaults');
  }

 
  #Slice methods return over hanging feature, so we could do this externally
  #If we are not returning over hanging features we should really handle whole chromosome
  #As we will need persistant data to handle over hangs
  #Do as whole chromosome anyway and chunk the Slices into windows which are divisible by all windows sizes
  #Nope, too much memmory, we either need to store here or generate the slice chunks in the caller and handle storing there
  #This is an adaptor after all let's do it here


  #Test for Rset here and whether we have status?
  #Currently just result_feature_set.
  #Need to implement roll back method in Helper

  #fetch this first to perform validation on slice and ResultSet??
  #my $rfs = $this->fetch_all_by_Slice_ResultSet($slice, $rset)
  #maybe want to tst if RESULT_FEATURE_SET before we do the query?


  #Need to as coord_system here.
  #Could rollback by slice(for all window sizes)


  if($rset->has_status('RESULT_FEATURE_SET')){
	throw('ResultSet('.$rset->name.') already has precomputed ResultFeatures stored, please rollback ResultFeature first');
	#Retrieving ResultFeatures here would end up retrieving them from the result_feature table
	#which is where we want to store them
	#They are only retrived from the probe/result/probe_feature table if they do not have this status.
	#REMEMBER TO SET THIS IN THE CALLING SCRIPT!

  }
  
  #Calulate sensible slice length based on window sizes
  my @wsizes = sort {$a <=> $b} @$window_sizes;
  
  #Current limit is 500KB, so let's start there
  #we should really start with a nice round number
  #so 


  my $multiplier = int(500000/$wsizes[$#wsizes]);
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
	
	warn "Could not find workable slice length for all window sizes, attempting to subset";
	
	#my %sorted_windows
	
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

	warn "largest chunk $largest_chunk with size $largest_size contains windows @largest_windows";

	my %remaining_windows = map {$_ => {}} @wsizes;
	delete $remaining_windows{'0'};#get rid of natural resolution as this will always work
	map { delete $remaining_windows{$_} } @largest_windows;

	#warn "Remaning windows are:\n".Data::Dumper::Dumper(\%remaining_windows);
	
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

	  warn "Found next chunk length $next_chunk contain remaining windows:\t@rwindows";

	#Now we want to cycle through all the set lengths which could contain the ones not in the first
	#so we need to
	}
	else{
	  
	  throw('Could not find workable slice length for remaining windows: '.
			join(', ', @rwindows));

	}

	warn "Need to write iterative sub for set definition";

  }
  else{
	warn "Found workable slice length $chunk_length for all window sizes";
	@{$chunk_windows{$chunk_length}} = keys(%workable_chunks);
	
  }

  #We need to subset the windows into workable sets
  #so 



  #Not lightweight as we will be storing them

  # Temporarily set the collection to be lightweight??? #Why?
  #my $old_value = $this->_lightweight();

  #if   ( defined($lightweight) ) { $this->_lightweight($lightweight) }
  #else                           { $this->_lightweight(1) }

  my (%counts, $store_natural);
  $store_natural = 1 if(grep/0/, @$window_sizes);

  
  foreach my $chunk_length(keys %chunk_windows){

	warn "Processing windows ".join(', ', @{$chunk_windows{$chunk_length}})." with chunk length $chunk_length";

	#Now walk through slice using slice length chunks and build all windows in each chunk
	my $in_slice     = 1;
	my $start_adj    = 0;
	my $end_adj      = $slice->start + $chunk_length - 1;
	my $region       = $slice->coord_system_name;
	my $seq_region_name  = $slice->seq_region_name;
	my $strand       = $slice->strand;
	my $rset_id      = $rset->dbID;
	
	my ($start, $end, $chunk_slice, $features, $bins);

	while($in_slice){

	  $start = $slice->start + $start_adj;
	  $end   = $slice->start + $end_adj;
	  
	  #Last chunk might not be the correct window length
	  #Hence why we should do this on whole chromosomes
	  #Force this or just warn, check seq_region slice
	  
	  if($end >= $slice->end){
		$end = $slice->end;
		$in_slice = 0;
	  }
	  

	  $chunk_slice = $slice_adaptor->fetch_by_region($region, $seq_region_name, $start, $end, $strand);
	  #warn "Processing slice ".$chunk_slice->name;

	  $features = $this->fetch_all_by_Slice_ResultSet($chunk_slice, $rset);

	  #warn "Got ".scalar(@$features);

	  #Shift chunk coords
	  if($in_slice){
		$start_adj += $chunk_length;
		$end_adj   += $chunk_length;
	  }

	
	  next if scalar(@$features) == 0;

	  #warn "\nFound ".scalar(@$features).' features';

	  
	  #This should return a hash of window size => bin array pairs
	  $bins = $this->_bin_features_by_window_sizes
				  ( 
				   -slice  => $slice,
				   -window_sizes  => $chunk_windows{$chunk_length},
				   -method => $method,
				   -features =>
				   $features,
				  );




	  #We need to handle strandedness of slice!?	  
	  my ($chunk_start, $chunk_end, $bin_score);
	  
	  

	  foreach my $wsize(keys %{$bins}){
		#warn "Got ".scalar(@{$bins->{$wsize}})." bins for window size $wsize";

		#We need to place feature back on original full length slice to store
		my $bin_start = $chunk_slice->start;
		my $bin_end   = $chunk_slice->start;
		$counts{$wsize} ||= 0;
	

		
	
		foreach my $bin_index(0..$#{$bins->{$wsize}}){
		  $bin_score = $bins->{$wsize}->[$bin_index];

		  
		  #next if ! $bin_score;#No we're no inc'ing the start ends for bins with no scores

		  $bin_end += $wsize;


		  if($bin_score){
			$counts{$wsize}++;	
			  
			#This is a little backwards as we are generating the object to store it
			#If we are aiming for speed the maybe we could also commodotise the store method
			#store by args or hash? store_fast?
			#Speed not essential for storing
			
			#Note: list ref passed
			
			#warn "storing $bin_start, $bin_end, $strand, $bin_score, undef, $rset->dbID, $wsize, $slice";
			
			$this->store(Bio::EnsEMBL::Funcgen::ResultFeature->new_fast
						 (
		  			  $bin_start, $bin_end, $strand, $bin_score, undef,#absent probe info
						  $rset->dbID, $wsize, $slice
						 ));
		  }
		  
		  $bin_start += $wsize;
		}
	  }



	  #Store all normal features in result_feature
	  if($store_natural){
		#warn "Storing natural resolution";

		foreach my $feature(@$features){
		  $counts{0}++;
		  
		  #warn "storing ".join(', ',	($feature->start, $feature->end, $feature->strand, $feature->score, 'undef', $rset->dbID, 0, $slice));


		  $this->store(Bio::EnsEMBL::Funcgen::ResultFeature->new_fast
					   (
						$feature->start, $feature->end, $feature->strand, $feature->score, undef,#absent probe info
						$rset->dbID, 0, $slice
						)); 
		}	
	  }
	  
	  
	  
	}



	#Turn off storing of natural resolution for next chunk length sets
	$store_natural = 0;

  }
  

  

  #print some counts here


  foreach my $wsize(keys %counts){
	warn "Stored ".$counts{$wsize}." for window size $wsize for ".$slice->name."\n";
  }

  #Return this counts hash so we can print/log from the caller, hence we don't print in here.
  
  return;
}


#This is a slight rewrite to calculate bins for a list of all window sizes
#Should only be called by store_window_bins_by_Slice_ResultSet
#This could be put in Collection and _bin_feature could call this to maintain interface

sub _bin_features_by_window_sizes{
  my $this = shift;
  my ( $slice, $window_sizes, $method, $features ) =
    rearrange( [ 'SLICE', 'WINDOW_SIZES', 'METHOD', 'FEATURES' ], @_ );


  if ( !defined($features) || !@{$features} ) { return {} }

  #No need to validate window sizes as done in the store_window_bins?
  #Do here anyway?

  #warn "Processing window sizes ".join(', ',@$window_sizes);


  #Set default to ResultFeature implementation?
  $method ||= 'average_score';

  #Should sub this to allow different implementations to use it
  $method = $this->validate_bin_method($method);


  my $slice_start = $slice->start();

  #my $bin_length = ( $slice->end() - $slice_start + 1 )/$nbins;
  #my @bins;

  #Set up some hashes to store data by window_size
  my (%bins, %nbins, %bin_counts);

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
  foreach my $wsize (@$window_sizes) {
	$nbins{$wsize}         = int($slice->length / $wsize);#should this be -1 as this will be an index?
	#no as int always rounds down
	#So nbins is actually the index of the bin not the 'number'

	#@{$bins{$wsize}}       = map( $first_bin, 1 .. $nbins{$wsize});
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
    
  #This failes for slices which are smaller than the chunk length;
  #Just do a little sanity test here to make sure the bin_lengths == wsizes
  #Then remove
  #foreach my $wsize (keys %bin_lengths) {
#	warn "wsize is $wsize with bin_length ".$bin_lengths{$wsize};

#	throw('Oh no, these should be the same') if($wsize != $bin_lengths{$wsize});

#  }
 
  my $feature_index = 0;

  #What is this?
  my @bin_masks;
  #my %start_bins;
  #my %end_bins;

  #Also store local starts and ends as we currently calc this several times
  my $lfeature_start;
  my $lfeature_end;


  #To remove skew from bins we need to clip the start and end 
  #if they do not have overlapping features.
  #So we need to store features in each bin and post process.
  #This would be easiest, but 2nd loop would slow it down.
  #Not so bad for storage, but let's try and keep it quick
  #So in line bin clipping requires awareness of when we have moved to a new bin
  #between features.  This holds for overlapping and non-overlapping features.
  #Once we have observed a gap we need to clip the end of the last bin and clip the start of the new bin.
  #This requires knowing the greatest end values from the last bin's feature.
  #what if two overlapping features had the same start and different end, would we see the longest last?
  #Check default slice_fetch sort


  foreach my $feature ( @{$features} ) {
	#Set up the bins for each window size

	foreach my $wsize (@$window_sizes) {
	
	  #We have already highjacked the object creation by here
	  #This is done in core BaseFeatureAdaptor
	  #We probably don't want to do this for ResultFeatures as we don't use the
	  #standard feature implementation
	  #we already use an array and we don't store the slice
	  #as this is already known by the caller
	  #and we always build on top level so we don't need to remap

	  #We do however need the slice to store?
	  #Do we?
	  #Yes, as we only store local starts when generating
	  #We need a store by Slice method?
	  #This will remove the need to inherit from Feature.
	  #These will need to be regenerated everytime we import a new build
	  #As we do with the probe_features themselves
	  #This also mean the result_feature status has to be associated with a coord_system_id
	  




	  #Which bins do the start and end lie in for this feature?
	  #Already dealing with local starts, so no slice subtraction


	  my $start_bin =  int(($feature->start ) / $wsize);
	  my $end_bin   =  int(($feature->end) / $wsize );
	
	
	  #my $start_bin =
	  #	int( ( $feature->[FEATURE_START] - $slice_start )/$bin_length );
	  
	  #  my $end_bin =
	  #	int( ( $feature->[FEATURE_END] - $slice_start )/$bin_length );
	  

	  #This might happen for last feature
	  #if (   $end_bins{$wsize} >= 	$nbins{$wsize} ) {
	  #	$end_bins{$wsize} =  $nbins{$wsize};

	  # }

	  $end_bin = $nbins{$wsize} if $end_bin > $nbins{$wsize};



	  #Slightly obfuscated code here
	  #Altho this should speed up generation
	  #by avoiding string comparisons.
	
	  #We should do default count processing for all methods as this is required for all no?
	  #No, just weight, coverage and average_score
	
=pod
	
	  if ( $method == 0 ) {
		throw('Not implemented for method for count/density');
		# ----------------------------------------------------------------
		# For 'count' and 'density'.
	  
		for ( my $bin_index = $start_bin ;
			  $bin_index <= $end_bin ;
			  ++$bin_index ) {
		  ++$bins[$bin_index];
		}

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

	  } elsif ( $method == 4 ) {
		# ----------------------------------------------------------------
		# For 'coverage'.
		
		#What exactly is this doing?
		#This is coverage of bin
		#Rather than coverage of feature as in fractional_count
	  
		
		my $feature_start = $feature->[FEATURE_START] - $slice_start;
		my $feature_end   = $feature->[FEATURE_END] - $slice_start;
		
		if ( !defined( $bin_masks[$start_bin] )
			 || ( defined( $bin_masks[$start_bin] )
				  && $bin_masks[$start_bin] != 1 ) ) {
		  # Mask the $start_bin from the start of the feature to the end
		  # of the bin, or to the end of the feature (whichever occurs
		  # first).
		  my $bin_start = int( $start_bin*$bin_length );
		  my $bin_end = int( ( $start_bin + 1 )*$bin_length - 1 );
		  for ( my $pos = $feature_start ;
				$pos <= $bin_end && $pos <= $feature_end ;
				++$pos ) {
			$bin_masks[$start_bin][ $pos - $bin_start ] = 1;
		  }
		}
		
		for ( my $bin_index = $start_bin + 1 ;
			  $bin_index <= $end_bin - 1 ;
			  ++$bin_index ) {
		  # Mark the middle bins between $start_bin and $end_bin as fully
		  # masked out.
		  $bin_masks[$bin_index] = 1;
		}
		
		if ( $end_bin != $start_bin ) {
		
		  if ( !defined( $bin_masks[$end_bin] )
			   || ( defined( $bin_masks[$end_bin] )
					&& $bin_masks[$end_bin] != 1 ) ) {
			# Mask the $end_bin from the start of the bin to the end of
			# the feature, or to the end of the bin (whichever occurs
			# first).
			my $bin_start = int( $end_bin*$bin_length );
			my $bin_end = int( ( $end_bin + 1 )*$bin_length - 1 );
			for ( my $pos = $bin_start ;
				  $pos <= $feature_end && $pos <= $bin_end ;
				  ++$pos ) {
			  $bin_masks[$end_bin][ $pos - $bin_start ] = 1;
			}
		  }
		  
		}
		
	  }							## end elsif ( $method == 4 )

=cut

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
	  } else {
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
  else{
	throw('Only accomodates average_score method');
  }


  #Could return bin_counts too summary reporting in zmenu
  #Could also do counting of specific type

  return \%bins;
} ## end sub _bin_features



#separated to allow addition of non-standard methods
#Could potentially add these in new
#and put this back in _bin_features


sub validate_bin_method{
  my ($self, $method) = @_;

  #Add average_score to avoid changing Collection.pm
  my $class = ref($self);
  ${$class::VALID_BINNING_METHODS}{'average_score'} = 5;



  
  #warn "Still can't access VALID_BINNING_METHODS";

  #foreach my $method_name(keys %{$class::VALID_BINNING_METHODS}){
#	warn "valid method is $method name";
#  }


  if ( ! exists( ${$class::VALID_BINNING_METHODS}{$method} ) ) {
    throw(
		  sprintf(
				  "Invalid binning method '%s', valid methods are:\n\t%s\n",
				  $method,
				  join( "\n\t", sort( keys(%{$self::VALID_BINNING_METHODS}) ) ) ) );
  }
  else{
	#warn "found valid method $method with index ".${$class::VALID_BINNING_METHODS}{$method};
  }
  
  return ${$class::VALID_BINNING_METHODS}{$method};
}



1;
