#
# EnsEMBL module for Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor
#
# Copyright (c) 2006 Ensembl
#
# You may distribute this module under the same terms as perl itself

=head1 NAME

Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor - An Base class for all
Funcgen FeatureAdaptors, redefines some methods to use the Funcgen DB

=head1 SYNOPSIS

Abstract class - should not be instantiated.  Implementation of
abstract methods must be performed by subclasses.

=head1 DESCRIPTION

This is a base adaptor for Funcgen feature adaptors. This base class is simply a way
of eliminating code duplication through the implementation of methods
common to all Funcgen feature adaptors.

=head1 CONTACT

Contact Ensembl development list for info: <ensembl-dev@ebi.ac.uk>

=head1 METHODS

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor;
use vars qw(@ISA @EXPORT);
use strict;


use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;
#use Bio::EnsEMBL::DBSQL::BaseAdaptor;#?
use Bio::EnsEMBL::Utils::Cache;
use Bio::EnsEMBL::Utils::Exception qw(warning throw deprecate stack_trace_dump);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);

@EXPORT = (@{$DBI::EXPORT_TAGS{'sql_types'}});

our $SLICE_FEATURE_CACHE_SIZE = 4;
our $MAX_SPLIT_QUERY_SEQ_REGIONS = 3;





#Added schema_build to feature cache key

=head2 fetch_all_by_Slice_constraint

  Arg [1]    : Bio::EnsEMBL::Slice $slice
               the slice from which to obtain features
  Arg [2]    : (optional) string $constraint
               An SQL query constraint (i.e. part of the WHERE clause)
  Arg [3]    : (optional) string $logic_name
               the logic name of the type of features to obtain
  Example    : $fs = $a->fetch_all_by_Slice_constraint($slc, 'perc_ident > 5');
  Description: Returns a listref of features created from the database which 
               are on the Slice defined by $slice and fulfill the SQL 
               constraint defined by $constraint. If logic name is defined, 
               only features with an analysis of type $logic_name will be 
               returned. 
  Returntype : listref of Bio::EnsEMBL::SeqFeatures in Slice coordinates
  Exceptions : thrown if $slice is not defined
  Caller     : Bio::EnsEMBL::Slice
  Status     : Stable

=cut

sub fetch_all_by_Slice_constraint {
  my($self, $slice, $constraint, $logic_name) = @_;

  my @result;

  if(!ref($slice) || !$slice->isa("Bio::EnsEMBL::Slice")) {
    throw("Bio::EnsEMBL::Slice argument expected.");
  }

  $constraint ||= '';
  $constraint = $self->_logic_name_to_constraint($constraint, $logic_name);

  #if the logic name was invalid, undef was returned
  return [] if(!defined($constraint));

  #check the cache and return if we have already done this query

  #Added schema_build to feature cache for EFG
  my $key = uc(join(':', $slice->name, $constraint, $self->db->_get_schema_build($slice->adaptor())));


  if(exists($self->{'_slice_feature_cache'}->{$key})) {
    return $self->{'_slice_feature_cache'}->{$key};
  }

  my $sa = $slice->adaptor();

  # Hap/PAR support: retrieve normalized 'non-symlinked' slices
  my @proj = @{$sa->fetch_normalized_slice_projection($slice)};

  if(@proj == 0) {
    throw('Could not retrieve normalized Slices. Database contains ' .
          'incorrect assembly_exception information.');
  }

  # Want to get features on the FULL original slice
  # as well as any symlinked slices

  # Filter out partial slices from projection that are on
  # same seq_region as original slice

  my $sr_id = $slice->get_seq_region_id();

  @proj = grep { $_->to_Slice->get_seq_region_id() != $sr_id } @proj;

  my $segment = bless([1,$slice->length(),$slice ],
                      'Bio::EnsEMBL::ProjectionSegment');
  push( @proj, $segment );


  # construct list of Hap/PAR boundaries for entire seq region
  my @bounds;
  my $ent_slice = $sa->fetch_by_seq_region_id($sr_id);
  $ent_slice    = $ent_slice->invert() if($slice->strand == -1);
  my @ent_proj  = @{$sa->fetch_normalized_slice_projection($ent_slice)};

  shift @ent_proj; # skip first
  @bounds = map {$_->from_start - $slice->start() + 1} @ent_proj;


  # fetch features for the primary slice AND all symlinked slices
  foreach my $seg (@proj) {
    my $offset = $seg->from_start();
    my $seg_slice  = $seg->to_Slice();

    my $features = $self->_slice_fetch($seg_slice, $constraint); ## NO RESULTS

    # if this was a symlinked slice offset the feature coordinates as needed
    if($seg_slice->name() ne $slice->name()) {

    FEATURE:
      foreach my $f (@$features) {
        if($offset != 1) {

          $f->{'start'} += $offset-1;
          $f->{'end'}   += $offset-1;
        }

        # discard boundary crossing features from symlinked regions
        foreach my $bound (@bounds) {
          if($f->{'start'} < $bound && $f->{'end'} >= $bound) {
            next FEATURE;
          }
        }

        $f->{'slice'} = $slice;
        push @result, $f;
      }
    }
    else {
      push @result, @$features;
    }
  }

  $self->{'_slice_feature_cache'}->{$key} = \@result;


  return \@result;
}


=head2 fetch_all_by_logic_name

  Arg [3]    : string $logic_name
               the logic name of the type of features to obtain
  Example    : $fs = $a->fetch_all_by_logic_name('foobar');
  Description: Returns a listref of features created from the database. 
               only features with an analysis of type $logic_name will be 
               returned. 
  Returntype : listref of Bio::EnsEMBL::SeqFeatures
  Exceptions : thrown if $logic_name
  Caller     : General
  Status     : Experimental

=cut

sub fetch_all_by_logic_name {
  my $self = shift;


  throw("Not implemented, is this relevant?");

  my $logic_name = shift || throw( "Need a logic_name" );
  my $constraint = $self->_logic_name_to_constraint( '',$logic_name );
  return $self->fetch_all( $constraint )
}

#
# helper function used by fetch_all_by_Slice_constraint method
#

#or maybe this is the one to change?

sub _slice_fetch {
	my $self = shift;
	my $slice = shift;
	my $orig_constraint = shift;
	
	my $slice_start  = $slice->start();
	my $slice_end    = $slice->end();
	my $slice_strand = $slice->strand();
	my $slice_cs     = $slice->coord_system();
	my $slice_seq_region = $slice->seq_region_name();
	
	#Need to get schema/data.version here from meta
	
	
	
	#get the synonym and name of the primary_table
	my @tabs = $self->_tables;
	my ($tab_name, $tab_syn) = @{$tabs[0]};
	
	#find out what coordinate systems the features are in
	my $mcc = $self->db->get_MetaCoordContainer();
	
	#reimplement and restrict to schema_build????????????????????????????????
	#This will select for all schema_builds
	#Should really restrict here, do below to avoid copying another module
	#or should we really do it in the MCC to enable other usage?
	my @feat_css = @{$mcc->fetch_all_CoordSystems_by_feature_type($tab_name)};
	
	my $asma = $self->db->get_AssemblyMapperAdaptor();
	my @features;
	#warn "jere are $self, @feat_css\n";

	# fetch the features from each coordinate system they are stored in
  COORD_SYSTEM: foreach my $feat_cs (@feat_css) {
		my $mapper;
		my @coords;
		my @ids;
		
		if($feat_cs->equals($slice_cs)) {#this now checks schema_build
			# no mapping is required if this is the same coord system
			
			my $max_len = $self->_max_feature_length() ||
			  $mcc->fetch_max_length_by_CoordSystem_feature_type($feat_cs,$tab_name);
			#should need to change this as we should have identified the correct EFG cs
		
			my $constraint = $orig_constraint;
			
			my $sr_id;
			if( $slice->adaptor() ) {
				$sr_id = $slice->adaptor()->get_seq_region_id($slice);
			} else {
				$sr_id = $self->db()->get_SliceAdaptor()->get_seq_region_id($slice);
			}
	

			

		
			$constraint .= " AND " if($constraint);
			$constraint .=
			  "${tab_syn}.seq_region_id = $sr_id AND " .
				"${tab_syn}.seq_region_start <= $slice_end AND " .
				  "${tab_syn}.seq_region_end >= $slice_start";
			

			#ensure correct coord_sys mapped to schema_build in $feat_cs->equals above
			$constraint .= " AND ${tab_syn}.coord_system_id = ".$feat_cs->dbID();


			if($max_len) {
				my $min_start = $slice_start - $max_len;
				$constraint .=
				  " AND ${tab_syn}.seq_region_start >= $min_start";
			}
			
  
	       
			my $fs = $self->generic_fetch($constraint,undef,$slice);
						
			# features may still have to have coordinates made relative to slice
			# start
			$fs = _remap($fs, $mapper, $slice);
			
			push @features, @$fs;
		} else {
			
			warn("Not yet implemented mapper for core to EFG DB");
			next;
			
			$mapper = $asma->fetch_by_CoordSystems($slice_cs, $feat_cs);
			
			next unless defined $mapper;
			
			# Get list of coordinates and corresponding internal ids for
			# regions the slice spans
			@coords = $mapper->map($slice_seq_region, $slice_start, $slice_end,
								   $slice_strand, $slice_cs);
			
			@coords = grep {!$_->isa('Bio::EnsEMBL::Mapper::Gap')} @coords;
			
			next COORD_SYSTEM if(!@coords);
			
			@ids = map {$_->id()} @coords;
			@ids = @{$asma->seq_regions_to_ids($feat_cs, \@ids)};
			
			# When regions are large and only partially spanned by slice
			# it is faster to to limit the query with start and end constraints.
			# Take simple approach: use regional constraints if there are less
			# than a specific number of regions covered.
			
			if(@coords > $MAX_SPLIT_QUERY_SEQ_REGIONS) {
				my $constraint = $orig_constraint;
				my $id_str = join(',', @ids);
				$constraint .= " AND " if($constraint);
				$constraint .= "${tab_syn}.seq_region_id IN ($id_str)";
				
				my $fs = $self->generic_fetch($constraint, $mapper, $slice);
			
				$fs = _remap($fs, $mapper, $slice);
			
				push @features, @$fs;

			} else {
				# do multiple split queries using start / end constraints
				
				my $max_len = $self->_max_feature_length() ||
				  $mcc->fetch_max_length_by_CoordSystem_feature_type($feat_cs,
																	 $tab_name);
				my $len = @coords;
				for(my $i = 0; $i < $len; $i++) {
					my $constraint = $orig_constraint;
					$constraint .= " AND " if($constraint);
					$constraint .=
					  "${tab_syn}.seq_region_id = "     . $ids[$i] . " AND " .
						"${tab_syn}.seq_region_start <= " . $coords[$i]->end() . " AND ".
						  "${tab_syn}.seq_region_end >= "   . $coords[$i]->start();
					
					if($max_len) {
						my $min_start = $coords[$i]->start() - $max_len;
						$constraint .=
						  " AND ${tab_syn}.seq_region_start >= $min_start";
					}
					
					my $fs = $self->generic_fetch($constraint,$mapper,$slice);
				
					$fs = _remap($fs, $mapper, $slice);
				
					push @features, @$fs;
				}
			}
		}
	} #COORD system loop

	return \@features;
}


####Redefined for Funcgen to trust the slice the feature is generated on
####and store on a non-core/dna i.e. slice enabled DB
####Also handles multiple coord systems, generating new coord_system_ids as appropriate
####Need to think about retrieval, i.e. mapping coord_system_id back to the core/dnadb
####version may not mean same seq_region_ids?  So need to use name somehow, name may not be unique!!
####Using the schema version may guarantee stability of seq_region_ids

#
# Helper function containing some common feature storing functionality
#
# Given a Feature this will return a copy (or the same feature if no changes 
# to the feature are needed) of the feature which is relative to the start
# of the seq_region it is on. The seq_region_id of the seq_region it is on
# is also returned.
#
# This method will also ensure that the database knows which coordinate
# systems that this feature is stored in.
#

sub _pre_store {
  my $self    = shift;
  my $feature = shift;

  if(!ref($feature) || !$feature->isa('Bio::EnsEMBL::Feature')) {
    throw('Expected Feature argument.');
  }

  $self->_check_start_end_strand($feature->start(),$feature->end(),
                                 $feature->strand());


  my $db = $self->db();
  my $slice = $feature->slice();

  if(!ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw('Feature must be attached to Slice to be stored.');
  }

  # make sure feature coords are relative to start of entire seq_region
  if($slice->start != 1 || $slice->strand != 1) {
	  throw("You must generate your feature on a slice starting at 1 with strand 1");
	  #move feature onto a slice of the entire seq_region
	  #$slice = $slice_adaptor->fetch_by_region($slice->coord_system->name(),
	  #                                         $slice->seq_region_name(),
	  #                                         undef, #start
	  #                                         undef, #end
	  #                                         undef, #strand
	  #                                         $slice->coord_system->version());
	  #$feature = $feature->transfer($slice);
	  #if(!$feature) {
	  #  throw('Could not transfer Feature to slice of ' .
	  #        'entire seq_region prior to storing');
	  #}
  }
  
  # Ensure this type of feature is known to be stored in this coord system.
  my $cs = $slice->coord_system;#from core/dnadb
 
  #retrieve corresponding Funcgen coord_system and set id in feature
  my $csa = $self->db->get_FGCoordSystemAdaptor();#had to call it FG as we were getting the core adaptor
  my $fg_cs = $csa->validate_coord_system($cs);
  $fg_cs = $csa->fetch_by_name_schema_build_version($cs->name(), $self->db->_get_schema_build($slice->adaptor()), $cs->version());
  $feature->coord_system_id($fg_cs->dbID());

  my ($tab) = $self->_tables();
  my $tabname = $tab->[0];


  #Need to do this for Funcgen DB
  my $mcc = $db->get_MetaCoordContainer();
  $mcc->add_feature_type($fg_cs, $tabname, $feature->length);

 # my $seq_region_id = $slice_adaptor->get_seq_region_id($slice);
  my $seq_region_id = $slice->get_seq_region_id();

  #would never get called as we're not validating against a different core DB
  #if(!$seq_region_id) {
  #  throw('Feature is associated with seq_region which is not in this dnadb.');
  #}

  return ($feature, $seq_region_id);
}



#
# helper function used to validate start/end/strand and 
# hstart/hend/hstrand etc.
#
sub _check_start_end_strand {
  my $self = shift;
  my $start = shift;
  my $end   = shift;
  my $strand = shift;

  #
  # Make sure that the start, end, strand are valid
  #
  if(int($start) != $start) {
    throw("Invalid Feature start [$start].  Must be integer.");
  }
  if(int($end) != $end) {
    throw("Invalid Feature end [$end]. Must be integer.");
  }
  if(int($strand) != $strand || $strand < -1 || $strand > 1) {
    throw("Invalid Feature strand [$strand]. Must be -1, 0 or 1.");
  }
  if($end < $start) {
    throw("Invalid Feature start/end [$start/$end]. Start must be less " .
          "than or equal to end.");
  }

  return 1;
}


#
# Given a list of features checks if they are in the correct coord system
# by looking at the first features slice.  If they are not then they are
# converted and placed on the slice.
#
sub _remap {
  my ($features, $mapper, $slice) = @_;

  #check if any remapping is actually needed

 
  if(@$features && (!$features->[0]->isa('Bio::EnsEMBL::Funcgen::Feature') ||
                    $features->[0]->slice == $slice)) {


    return $features;
  }

      


  throw("Not yet fully implemented _remap for EFG");

  #remapping has not been done, we have to do our own conversion from
  #to slice coords

  my @out;

  my $slice_start = $slice->start();
  my $slice_end   = $slice->end();
  my $slice_strand = $slice->strand();
  my $slice_cs    = $slice->coord_system();

  my ($seq_region, $start, $end, $strand);

  my $slice_seq_region = $slice->seq_region_name();

  foreach my $f (@$features) {
    #since feats were obtained in contig coords, attached seq is a contig
    my $fslice = $f->slice();
    if(!$fslice) {
      throw("Feature does not have attached slice.\n");
    }
    my $fseq_region = $fslice->seq_region_name();
    my $fcs = $fslice->coord_system();

	#This should never be the case for efg??
	#Note:  Using non EFG equals sub here, so no schema_build check
    if(!$slice_cs->equals($fcs)) {
      #slice of feature in different coord system, mapping required

      ($seq_region, $start, $end, $strand) =
        $mapper->fastmap($fseq_region,$f->start(),$f->end(),$f->strand(),$fcs);

      # undefined start means gap
      next if(!defined $start);
    } else {
      $start      = $f->start();
      $end        = $f->end();
      $strand     = $f->strand();
      $seq_region = $f->slice->seq_region_name();
    }

    # maps to region outside desired area
    next if ($start > $slice_end) || ($end < $slice_start) || 
      ($slice_seq_region ne $seq_region);

    #shift the feature start, end and strand in one call
    if($slice_strand == -1) {
      $f->move( $slice_end - $end + 1, $slice_end - $start + 1, $strand * -1 );
    } else {
      $f->move( $start - $slice_start + 1, $end - $slice_start + 1, $strand );
    }

    $f->slice($slice);

    push @out,$f;
  }

  return \@out;
}


#
# Given a logic name and an existing constraint this will
# add an analysis table constraint to the feature.  Note that if no
# analysis_id exists in the columns of the primary table then no
# constraint is added at all
#
sub _logic_name_to_constraint {
  my $self = shift;
  my $constraint = shift;
  my $logic_name = shift;

  return $constraint if(!$logic_name);

  #make sure that an analysis_id exists in the primary table
  my ($prim_tab) = $self->_tables();
  my $prim_synonym = $prim_tab->[1];

  my $found_analysis=0;
  foreach my $col ($self->_columns) {
    my ($syn,$col_name) = split(/\./,$col);
    next if($syn ne $prim_synonym);
    if($col_name eq 'analysis_id') {
      $found_analysis = 1;
      last;
    }
  }

  if(!$found_analysis) {
    warning("This feature is not associated with an analysis.\n" .
            "Ignoring logic_name argument = [$logic_name].\n");
    return $constraint;
  }

  my $aa = $self->db->get_AnalysisAdaptor();
  my $an = $aa->fetch_by_logic_name($logic_name);

  if(!$an) {
    return undef;
  }

  my $an_id = $an->dbID();

  $constraint .= ' AND' if($constraint);
  $constraint .= " ${prim_synonym}.analysis_id = $an_id";
  return $constraint;
}


=head2 store

  Arg [1]    : list of Bio::EnsEMBL::SeqFeature
  Example    : $adaptor->store(@feats);
  Description: ABSTRACT  Subclasses are responsible for implementing this 
               method.  It should take a list of features and store them in 
               the database.
  Returntype : none
  Exceptions : thrown method is not implemented by subclass
  Caller     : general
  Status     : At Risk
             : throws if called.

=cut

sub store{
  my $self = @_;

  throw("Abstract method store not defined by implementing subclass\n");
}


=head2 remove

  Arg [1]    : A feature $feature 
  Example    : $feature_adaptor->remove($feature);
  Description: This removes a feature from the database.  The table the
               feature is removed from is defined by the abstract method
               _tablename, and the primary key of the table is assumed
               to be _tablename() . '_id'.  The feature argument must 
               be an object implementing the dbID method, and for the
               feature to be removed from the database a dbID value must
               be returned.
  Returntype : none
  Exceptions : thrown if $feature arg does not implement dbID(), or if
               $feature->dbID is not a true value
  Caller     : general
  Status     : Stable

=cut


sub remove {
  my ($self, $feature) = @_;

  throw("Not yet implemented for efg");

  if(!$feature || !ref($feature) || !$feature->isa('Bio::EnsEMBL::Feature')) {
    throw('Feature argument is required');
  }

  if(!$feature->is_stored($self->db)) {
    throw("This feature is not stored in this database");
  }

  my @tabs = $self->_tables;
  my ($table) = @{$tabs[0]};

  my $sth = $self->prepare("DELETE FROM $table WHERE ${table}_id = ?");
  $sth->bind_param(1,$feature->dbID,SQL_INTEGER);
  $sth->execute();

  #unset the feature dbID ad adaptor
  $feature->dbID(undef);
  $feature->adaptor(undef);

  return;
}


=head2 remove_by_Slice

  Arg [1]    : Bio::Ensembl::Slice $slice
  Example    : $feature_adaptor->remove_by_Slice($slice);
  Description: This removes features from the database which lie on a region
               represented by the passed in slice.  Only features which are
               fully contained by the slice are deleted; features which overlap
               the edge of the slice are not removed.
               The table the features are removed from is defined by
               the abstract method_tablename.
  Returntype : none
  Exceptions : thrown if no slice is supplied
  Caller     : general
  Status     : Stable

=cut

sub remove_by_Slice {
  my ($self, $slice) = @_;

  if(!$slice || !ref($slice) || !$slice->isa('Bio::EnsEMBL::Slice')) {
    throw("Slice argument is required");
  }

  throw("Not yet implemented for efg");


  my @tabs = $self->_tables;
  my ($table_name) = @{$tabs[0]};

  my $seq_region_id = $self->db->get_SliceAdaptor->get_seq_region_id($slice);
  my $start = $slice->start();
  my $end   = $slice->end();

  #
  # Delete only features fully on the slice, not overlapping ones
  #
  my $sth = $self->prepare("DELETE FROM $table_name " .
                           "WHERE seq_region_id = ? " .
                           "AND   seq_region_start >= ? " .
                           "AND   seq_region_end <= ?");

  $sth->bind_param(1,$seq_region_id,SQL_INTEGER);
  $sth->bind_param(2,$start,SQL_INTEGER);
  $sth->bind_param(3,$end,SQL_INTEGER);
  $sth->execute();
  $sth->finish();
}


#
# Internal function. Allows the max feature length which is normally
# retrieved from the meta_coord table to be overridden.  This allows
# for some significant optimizations to be put in when it is known
# that requested features will not be over a certain size.
#
sub _max_feature_length {
  my $self = shift;
  return $self->{'_max_feature_length'} = shift if(@_);
  return $self->{'_max_feature_length'};
}





1;


