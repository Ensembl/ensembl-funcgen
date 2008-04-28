# $Id: ResultFeature.pm,v 1.1 2008-04-28 19:12:12 nj1 Exp $

package Bio::EnsEMBL::Collection::Exon;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument  ('rearrange');
use Bio::EnsEMBL::Utils::Exception ('throw');

use base( 'Bio::EnsEMBL::Collection',
          'Bio::EnsEMBL::DBSQL::ExonAdaptor' );

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
    @columns[ 5 .. $#columns ] = map( 1, 5 .. $#columns );
  }

  return @columns;
}

sub _default_where_clause {
  my ($this) = @_;

  if ( $this->_lightweight() ) {
    return '';
  }

  return $this->SUPER::_default_where_clause();
}




sub store_window_bins_by_Slice_ResultSet {
  my $this = shift;
  my ( $slice, $rset, $window_sizes, $logic_name ) =
    rearrange( [ 'SLICE', 'RESULT_SET', 'WINDOW_SIZES', 'LOGIC_NAME' ], @_ );


  my $slice_adaptor = $this->db->dnadb->get_SliceAdaptor;
  my $method='coverage';#average
  #Hardcoded method to average.  What about those probe which cross boundaries?
  #Will this be returned?
  throw('Must supply at least one window size as an arrayref') if scalar(@$window_sizes) <= 0;
  
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


  if($rset->has_status('RESULT_FEATURE_SET')){
	throw('ResultSet('.$rset->name.') already has precomputed ResultFeatures stored, please rollback ResultFeature first');
	#Retrieving ResultFeatures here would end up retrieving them from the result_feature table
	#which is where we want to store them
	#They are only retrived from the probe/result/probe_feature table if they do not have this status.
	#REMEMBER TO SET THIS IN THE CALLING SCRIPT!

  }
  
  #Calulate sensible slice length based on window sizes
  my @wsizes = sort $a <=> $b, @$window_sizes;
  
  #Current limit is 500KB, so let's start there
  my ($multiplier) = split/\./, (500000/$wsizes[$#wizes]);
  my $chunk_length = $multiplier * $wsizes[$#wizes];
  warn "multiplier is $multiplier";
  my $not_divisible = 1;
  
  while($not_divisible && $slice_length != 0){
	$not_divisible = 0;
	
	foreach my $wsize(@wsizes){
	  next if $wsize == 0;#Special wsize for normal data
	  
	  #Set not divisible if modulus is true
	  $not_divisible = 1 if($slice_length % $wsize);
	}
	
	#Gradually shrink the length until we find a workable slice length for all windows
	$chunk_length -= $wsizes[$#wizes] if $not_divisible;
  }
  
  throw('Could not find workable slice length below 500kb for windows: '.
		join(', ', @wsizes)) if $slice_length == 0;
  
  warn "Found workable slice length $slice_length";


  #Not lightweight as we will be storing them

  # Temporarily set the collection to be lightweight??? #Why?
  #my $old_value = $this->_lightweight();

  #if   ( defined($lightweight) ) { $this->_lightweight($lightweight) }
  #else                           { $this->_lightweight(1) }


  #Now walk through slice using slice length chunks and build all windows in each chunk
  my $in_slice     = 1;
  my $start_adj    = 0;
  my $end_adj      = $slice->seq_region_start + $chunk_length - 1;
  my $region       = $slice->coord_system_name;
  my $region_name  = $slice->seq_region_name;
  my $strand       = $slice->strand;

  my ($start, $end, $chunk_slice, $bins);

  while($in_slice){

	$start = $slice->seq_region_start + $start_adj;
	$end   = $slice->seq_region_start + $end_adj;

	#Last chunk might not be the correct window length
	#Hence why we should do this on whole chromosomes
	#Force this or just warn, check seq_region slice

	if($end >= $slice->seq_region_end){
	  $end = $slice->seq_region_end;
	  $in_slice = 0;
	}

	$chunk_slice = $slice_adaptor->fetch_by_region($region, $seq_region_name, $start, $end, $strand);

	#This should return a hash of window size => bin array pairs
	$bins =
	  $this->_bin_features_by_window_size( -slice  => $slice,
										   -window_sizes  => $window_sizes,
										   -method => $method,
										   -features =>
										   $this->fetch_all_by_Slice_ResultSet($chunk_slice, $rset),
										 );
  
	#Now store

	foreach my $wsize(keys %{$bins}){
	  
	  foreach my $result_feat_bin(@{$bins->{$wsize}}){

		#Do some counting here
		
		#This is a little backwards as we are generating the object to store it
		#If we are aiming for speed the maybe we could also commodotise the store method
		#store by args or hash? store_fast?
		$self->store(Bio::EnsEMBL::Funcgen::ResultFeature->new(
															   #(seq_region_start, seq_region_end, score, result_set_id, seq_region_id, seq_region_strand, window_size
					 
															  )
					);
		
	
	  }

	}

	if($in_slice){
	  $start_adj += $chunk_length;
	  $end_adj   += $chunk_length
	}

  }

  #print some counts here

  return;

}

  # Reset the lightweightness to whatever it was before.
  #$this->_lightweight($old_value);

  return $bins;
}

#This is a slight rewrite to calculate bins for a list of all window sizes
#Should only be called by store_window_bins_by_Slice_ResultSet

sub _bin_features_by_window_size{
  my $this = shift;
  my ( $slice, $window_sizes, $method_name, $features ) =
    rearrange( [ 'SLICE', 'WINDOW_SIZES', 'METHOD', 'FEATURES' ], @_ );

  if ( !defined($features) || !@{$features} ) { return [] }

  #if ( !defined($nbins) ) {
  #  throw('Missing NBINS argument');
  #} elsif ( $nbins <= 0 ) {
  #  throw('Negative or zero NBINS argument');
  #}

  #No need to validate window sizes as done in the caller

  $method_name ||= 'count';
  if ( !exists( $VALID_BINNING_METHODS{$method_name} ) ) {
    throw(
           sprintf(
              "Invalid binning method '%s', valid methods are:\n\t%s\n",
              $method_name,
              join( "\n\t", sort( keys(%VALID_BINNING_METHODS) ) ) ) );
  }
  my $method = $VALID_BINNING_METHODS{$method_name};

  my $slice_start = $slice->start();

  #my $bin_length = ( $slice->end() - $slice_start + 1 )/$nbins;
  #my @bins;

  #Set up some hashes to store data by window_size
  my (%bin_lengths, %window_bins);

  if ( $method == 0 ||    # 'count' or 'density'
       $method == 3 ||    # 'fractional_count' or 'weight'
       $method == 4       # 'coverage'
    ){
    # For binning methods where each bin contain numerical values.
	$first_bin =  = 0;
  } 
  else {
	# For binning methods where each bin does not contain numerical
    # values.
	$first_bin = undef;
  }

  foreach my $wsize(@$window_sizes){
	my ($nbins) = split/\./,($slice->length / $wsize);
	$bin_lengths{$wsize} = ( $slice->end() - $slice_start + 1 )/$nbins;#This should be wsize!!!
	$window_bins{$wsize} = map( $first_bin, 1 .. $nbins );
  }
  

  #Just do a little sanity test here to make sure the bin_lengths == wsizes
  #Then remove

  foreach my $wsize(keys %bin_lengths){
	warn "wsize is $wsize with bin_length ".$bin_length{$wsize};

	throw('Oh no, these should be the same');

  }
  
  warn 'Can remove wsize versus bin_length chceck and change assignment calc';

  my $feature_index = 0;
  my @bin_masks;

  foreach my $feature ( @{$features} ) {
    my $start_bin =
      int( ( $feature->[FEATURE_START] - $slice_start )/$bin_length );
    my $end_bin =
      int( ( $feature->[FEATURE_END] - $slice_start )/$bin_length );

    if ( $end_bin >= $nbins ) {
      # This might happen for the very last entry.
      $end_bin = $nbins - 1;
    }

    if ( $method == 0 ) {
      # ----------------------------------------------------------------
      # For 'count' and 'density'.

      for ( my $bin_index = $start_bin ;
            $bin_index <= $end_bin ;
            ++$bin_index )
      {
        ++$bins[$bin_index];
      }

    } elsif ( $method == 1 ) {
      # ----------------------------------------------------------------
      # For 'indices' and 'index'

      for ( my $bin_index = $start_bin ;
            $bin_index <= $end_bin ;
            ++$bin_index )
      {
        push( @{ $bins[$bin_index] }, $feature_index );
      }

      ++$feature_index;

    } elsif ( $method == 2 ) {
      # ----------------------------------------------------------------
      # For 'features' and 'feature'.

      for ( my $bin_index = $start_bin ;
            $bin_index <= $end_bin ;
            ++$bin_index )
      {
        push( @{ $bins[$bin_index] }, $feature );
      }

    } elsif ( $method == 3 ) {
      # ----------------------------------------------------------------
      # For 'fractional_count' and 'weight'.

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
              ++$bin_index )
        {
          $bins[$bin_index] += $bin_length/$feature_length;
        }

        # The last bin...
        $bins[$end_bin] +=
          ( ( $feature->[FEATURE_END] - $slice_start ) -
            $end_bin*$bin_length +
            1 )/$feature_length;

      } ## end else [ if ( $start_bin == $end_bin)

    } elsif ( $method == 4 ) {
      # ----------------------------------------------------------------
      # For 'coverage'.

      my $feature_start = $feature->[FEATURE_START] - $slice_start;
      my $feature_end   = $feature->[FEATURE_END] - $slice_start;

      if ( !defined( $bin_masks[$start_bin] )
           || ( defined( $bin_masks[$start_bin] )
                && $bin_masks[$start_bin] != 1 ) )
      {
        # Mask the $start_bin from the start of the feature to the end
        # of the bin, or to the end of the feature (whichever occurs
        # first).
        my $bin_start = int( $start_bin*$bin_length );
        my $bin_end = int( ( $start_bin + 1 )*$bin_length - 1 );
        for ( my $pos = $feature_start ;
              $pos <= $bin_end && $pos <= $feature_end ;
              ++$pos )
        {
          $bin_masks[$start_bin][ $pos - $bin_start ] = 1;
        }
      }

      for ( my $bin_index = $start_bin + 1 ;
            $bin_index <= $end_bin - 1 ;
            ++$bin_index )
      {
        # Mark the middle bins between $start_bin and $end_bin as fully
        # masked out.
        $bin_masks[$bin_index] = 1;
      }

      if ( $end_bin != $start_bin ) {

        if ( !defined( $bin_masks[$end_bin] )
             || ( defined( $bin_masks[$end_bin] )
                  && $bin_masks[$end_bin] != 1 ) )
        {
          # Mask the $end_bin from the start of the bin to the end of
          # the feature, or to the end of the bin (whichever occurs
          # first).
          my $bin_start = int( $end_bin*$bin_length );
          my $bin_end = int( ( $end_bin + 1 )*$bin_length - 1 );
          for ( my $pos = $bin_start ;
                $pos <= $feature_end && $pos <= $bin_end ;
                ++$pos )
          {
            $bin_masks[$end_bin][ $pos - $bin_start ] = 1;
          }
        }

      }

    } ## end elsif ( $method == 4 )

  } ## end foreach my $feature ( @{$features...

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

  return \@bins;
} ## end sub _bin_features



1;
