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
use Bio::EnsEMBL::Utils::Cache;
use Bio::EnsEMBL::Utils::Exception qw(warning throw deprecate stack_trace_dump);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);

@EXPORT = (@{$DBI::EXPORT_TAGS{'sql_types'}});

our $SLICE_FEATURE_CACHE_SIZE = 4;
our $MAX_SPLIT_QUERY_SEQ_REGIONS = 3;


#need to wrap _generic_fetch to always generate the seq_region_cache other wise non-slice based
#fetch methods fail

=head2 generic_fetch

  Arg [1]    : (optional) string $constraint
               An SQL query constraint (i.e. part of the WHERE clause)
  Arg [2]    : (optional) Bio::EnsEMBL::AssemblyMapper $mapper
               A mapper object used to remap features
               as they are retrieved from the database
  Arg [3]    : (optional) Bio::EnsEMBL::Slice $slice
               A slice that features should be remapped to
  Example    : $fts = $a->generic_fetch('contig_id in (1234, 1235)', 'Swall');
  Description: Wrapper method for core BaseAdaptor, build seq_region cache for features
  Returntype : ARRAYREF of Bio::EnsEMBL::SeqFeature in contig coordinates
  Exceptions : none
  Caller     : FeatureAdaptor classes
  Status     : at risk

=cut

sub generic_fetch {
  my $self = shift;

  #, $constraint, $mapper, $slice) = @_;
  #build seq_region cache here once for entire query
  $self->build_seq_region_cache();

  return $self->SUPER::generic_fetch(@_);
}



#Added schema_buildo feature cache key
#Do not remove

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

  #add cs_id here to amke sure we're definitely getting
  #feature from the correct seq_region, as seq_region_ids can be re-used for
  #different seq_region between DB with different default assemblies

  $constraint ||= '';#need to addcs_id here and coord_system_id='.$slice;
  

  my $fg_cs = $self->db->get_FGCoordSystemAdaptor->fetch_by_name(
																$slice->coord_system->name(), 
																$slice->coord_system->version()
															   );


  if(! defined $fg_cs){
	warn "No CoordSystem present for ".$slice->coord_system->name().":".$slice->coord_system->version();
	return \@result;
  }

  #build seq_region cache here once for entire query
  $self->build_seq_region_cache($slice);#, $fg_cs);


  my @tables = $self->_tables;
  my (undef, $syn) = @{$tables[0]};

  #$constraint .= ' AND ' if $constraint;#constraint can be empty string
  #$constraint .= " ${syn}.coord_system_id=".$fg_cs->dbID();
  $constraint = $self->_logic_name_to_constraint($constraint, $logic_name);

  #if the logic name was invalid, undef was returned
  return [] if(!defined($constraint));

  #check the cache and return if we have already done this query

  #Added schema_build to feature cache for EFG
  my $key = uc(join(':', $slice->name, $constraint, $self->db->_get_schema_build($slice->adaptor->db())));


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

    my $features = $self->_slice_fetch($seg_slice, $constraint); ## NO RESULTS? This is a problem with the cs->equals method?


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


#do we even need to have the coord system?
#so lon as we are only using one schema build i.e. one dnadb defualt = current
#slice and fg_cs are optional
#need to look at this

sub build_seq_region_cache{
  my ($self, $slice) = @_;

  #if(defined $fg_cs && !(ref($fg_cs) && $fg_cs->isa('Bio::EnsEMBL::Funcgen::CoordSystem'))){
#	throw('Optional argument must be a Bio::EnsEMBL::Funcgen::CoordSystem');
#  }else{
#	$fg_cs = $self->db->get_FGCoordSystemAdaptor->fetch_by_name(
#																$slice->coord_system->name(), 
#																$slice->coord_system->version()
#															   );
#  }

  if(defined $slice){
	throw('Optional argument must be a Bio::EnsEMBL::Slice') if(! ( ref($slice) && $slice->isa('Bio::EnsEMBL::Slice')));
  }


  my $dnadb = (defined $slice) ? $slice->adaptor->db() : $self->db->dnadb();
  my $schema_build = $self->db->_get_schema_build($dnadb);

  #my $sql = 'SELECT core_seq_region_id, seq_region_id from seq_region where coord_system_id='.
#	$fg_cs->dbID().' and schema_build="'.$schema_build.'"';
  my $sql = 'SELECT core_seq_region_id, seq_region_id from seq_region where schema_build="'.$schema_build.'"';

  
  $self->{'seq_region_cache'} = {};
  $self->{'core_seq_region_cache'} = {};
  %{$self->{'seq_region_cache'}} = map @$_, @{$self->db->dbc->db_handle->selectall_arrayref($sql)};

  
  #this return nothing for a new schema_build!!
  #How can

  #now reverse cache
  
  foreach my $csr_id (keys %{$self->{'seq_region_cache'}}){

	$self->{'core_seq_region_cache'}->{$self->{'seq_region_cache'}->{$csr_id}} = $csr_id;
  }

  return;
}


sub get_seq_region_id_by_Slice{
  my ($self, $slice) = @_;

  if(! ($slice && ref($slice) && $slice->isa("Bio::EnsEMBL::Slice"))){
	throw('You must provide a valid Bio::EnsEMBL::Slice');
  }

  my $sr_id;

  if( $slice->adaptor() ) {
	$sr_id = $slice->adaptor()->get_seq_region_id($slice);
  } else {
	$sr_id = $self->db()->get_SliceAdaptor()->get_seq_region_id($slice);
  }

  #should we wanr or thro here if not exists?

  return (exists $self->{'seq_region_cache'}->{$sr_id}) ? $self->{'seq_region_cache'}->{$sr_id} : undef;
}

sub get_core_seq_region_id{
  my ($self, $fg_sr_id) = @_;

  return $self->{'core_seq_region_cache'}->{$fg_sr_id};
}

=head2 _pre_store

  Arg [1]    : Bio::EnsEMBL::Feature
  Example    : $fs = $a->fetch_all_by_Slice_constraint($slc, 'perc_ident > 5');
  Description: Helper function containing some common feature storing functionality
               Given a Feature this will return a copy (or the same feature if no changes 
	       to the feature are needed) of the feature which is relative to the start
               of the seq_region it is on. The seq_region_id of the seq_region it is on
               is also returned.  This method will also ensure that the database knows which coordinate
               systems that this feature is stored in.  This supercedes teh core method, to trust the
               slice the feature has been generated on i.e. from the dnadb.  Also handles multi-coordsys
               aspect, generating new coord_system_ids as appropriate
  Returntype : Bio::EnsEMBL::Feature and the seq_region_id it is mapped to
  Exceptions : thrown if $slice is not defined
  Caller     : Bio::EnsEMBL::"Type"FeatureAdaptors
  Status     : At risk

=cut


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
  my $fg_cs = $csa->validate_and_store_coord_system($cs);

  
  $fg_cs = $csa->fetch_by_name($cs->name(), $cs->version());

  #removed as we don't want this to slow down import
  #my $sbuild = $self->db->_get_schema_build($slice->adaptor->db());

  #if(! $fg_cs->contains_schema_build($sbuild)){

	#warn "Adding new schema build $sbuild to CoordSystem\n";

	#$fg_cs->add_core_coord_system_info(
		#							   -RANK                 => $cs->rank(), 
		#							   -SEQ_LVL              => $cs->is_sequence_level(), 
		#							   -DEFAULT_LVL          => $cs->is_default(), 
		#							   -SCHEMA_BUILD         => $sbuild, 
		#							   -CORE_COORD_SYSTEM_ID => $cs->dbID(),
		#							   -IS_STORED            => 0,
		#							  );
	#$fg_cs = $csa->store($fg_cs);
  #}

  #$feature->coord_system_id($fg_cs->dbID());

  my ($tab) = $self->_tables();
  my $tabname = $tab->[0];


  #Need to do this for Funcgen DB
  my $mcc = $db->get_MetaCoordContainer();
  $mcc->add_feature_type($fg_cs, $tabname, $feature->length);


  #build seq_region cache here once for entire query
  $self->build_seq_region_cache($slice);#, $fg_cs);

  #Now need to check whether seq_region is already stored
  my $seq_region_id = $self->get_seq_region_id_by_Slice($slice);

  if(! $seq_region_id){
	#check whether we have an equivalent seq_region_id

	my $sql = 'SELECT seq_region_id from seq_region where coord_system_id='.$fg_cs->dbID().
	  ' and name="'.$slice->seq_region_name.'"';

  	($seq_region_id) = $self->db->dbc->db_handle->selectrow_array($sql);
	my $schema_build = $self->db->_get_schema_build($slice->adaptor->db());

	#No compararble seq_region
	if(! $seq_region_id){
	  $sql = 'INSERT into seq_region(name, coord_system_id, core_seq_region_id, schema_build) '.
		'values("'.$slice->seq_region_name().'", '.$fg_cs->dbID().', '.$slice->get_seq_region_id().', "'.$schema_build.'")';
	}else{#Add to comparable seq_region
	  $sql = 'INSERT into seq_region(seq_region_id, name, coord_system_id, core_seq_region_id, schema_build) '.
		'values('.$seq_region_id.', "'.$slice->seq_region_name().'", '.$fg_cs->dbID().', '.$slice->get_seq_region_id().', "'.$schema_build.'")';
	}

	my $sth = $self->prepare($sql);
	$sth->execute();
	$seq_region_id =  $sth->{'mysql_insertid'};
  }

  
  #my $seq_region_id = $slice->get_seq_region_id();

  #would never get called as we're not validating against a different core DB
  #if(!$seq_region_id) {
  #  throw('Feature is associated with seq_region which is not in this dnadb.');
  #}
  

  #why are returning seq_region id here..is thi snot just in the feature slice?
  #This is actually essential as the seq_region_ids are not stored in the slice
  #retrieved from slice adaptor
  return ($feature, $seq_region_id);
}


#This is causing problems with remapping

sub _slice_fetch {
  my $self = shift;
  my $slice = shift;
  my $orig_constraint = shift;

  my $slice_start  = $slice->start();
  my $slice_end    = $slice->end();
  my $slice_strand = $slice->strand();
  my $slice_cs     = $slice->coord_system();
  #my $slice_seq_region = $slice->seq_region_name();


  #### Here!!! We Need to translate the seq_regions IDs to efg seq_region_ids
  #we need to fetch the seq_region ID based on the coord_system id and the name
  #we don't want to poulate with the eFG seq_region_id, jsut the core one, as we need to maintain core info in the slice.
  
  #we need to cache the seq_region_id mappings for the coord_system and schema_build, 
  #so we don't do it for every feature in a slice
  #should we just reset it in temporarily in fetch above, other wise we lose cache between projections
  #should we add seq_region cache to coord_system/adaptor?
  #Then we can access it from here and never change the slice seq_region_id
  #how are we accessing the eFG CS?  As we only have the core CS from the slice here?
  #No point in having in eFG CS if we're just having to recreate it for everyquery
  #CSA already has all CS's cached, so query overhead would be mininal
  #But cache overhead for all CS seq_region_id caches could be quite large.
  #we're going to need it in the _pre_store, so having it in the CS would be natural there
  #can we just cache dynamically instead?
  #we're still calling for csa rather than accessing from cache
  #could we cache the csa in the base feature adaptor?
  

  

  #my $slice_seq_region_id = $slice->get_seq_region_id();

  #get the synonym and name of the primary_table
  my @tabs = $self->_tables;
  my ($tab_name, $tab_syn) = @{$tabs[0]};

  #find out what coordinate systems the features are in
  my $mcc = $self->db->get_MetaCoordContainer();
  my @feat_css=();


  my $mca = $self->db->get_MetaContainer();
  my $value_list = $mca->list_value_by_key( $tab_name."build.level" );
  if( @$value_list and $slice->is_toplevel()) {   
    push @feat_css, $slice_cs;
  }
  else{
    @feat_css = @{$mcc->fetch_all_CoordSystems_by_feature_type($tab_name)};
  }

  my $asma = $self->db->get_AssemblyMapperAdaptor();
  my @features;

  # fetch the features from each coordinate system they are stored in
 COORD_SYSTEM: foreach my $feat_cs (@feat_css) {

    my $mapper;
    my @coords;
    my @ids;

    if($feat_cs->equals($slice_cs)) {
	  #warn "got match";
      # no mapping is required if this is the same coord system

      my $max_len = $self->_max_feature_length() ||
        $mcc->fetch_max_length_by_CoordSystem_feature_type($feat_cs,$tab_name);

      my $constraint = $orig_constraint;

      #my $sr_id;
      #if( $slice->adaptor() ) {
	  #$sr_id = $self->get_seq_region_id_by_Slice($slice->adaptor()->get_seq_region_id($slice));
      #} else {
	  #	$sr_id = $self->get_seq_region_id_by_Slice($self->db()->get_SliceAdaptor()->get_seq_region_id($slice));
	  #  }

	  my $sr_id = $self->get_seq_region_id_by_Slice($slice);


      $constraint .= " AND " if($constraint);
      $constraint .=
          "${tab_syn}.seq_region_id = $sr_id AND " .
          "${tab_syn}.seq_region_start <= $slice_end AND " .
          "${tab_syn}.seq_region_end >= $slice_start";


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
    } 

	#can't do remapping yet as AssemblyMapper expects a core CS
	#change AssemblyMapper
	#or do we just create a core CS just for the remap and convert back when done

	#else {
    #  $mapper = $asma->fetch_by_CoordSystems($slice_cs, $feat_cs);

   #   next unless defined $mapper;

      # Get list of coordinates and corresponding internal ids for
      # regions the slice spans
   #   @coords = $mapper->map($slice_seq_region, $slice_start, $slice_end,
    #                         $slice_strand, $slice_cs);

    #  @coords = grep {!$_->isa('Bio::EnsEMBL::Mapper::Gap')} @coords;

    #  next COORD_SYSTEM if(!@coords);

     # @ids = map {$_->id()} @coords;
  ##coords are now id rather than name
  ##      @ids = @{$asma->seq_regions_to_ids($feat_cs, \@ids)};

      # When regions are large and only partially spanned by slice
      # it is faster to to limit the query with start and end constraints.
      # Take simple approach: use regional constraints if there are less
      # than a specific number of regions covered.

    #  if(@coords > $MAX_SPLIT_QUERY_SEQ_REGIONS) {
    #    my $constraint = $orig_constraint;
    #    my $id_str = join(',', @ids);
    #    $constraint .= " AND " if($constraint);
    #    $constraint .= "${tab_syn}.seq_region_id IN ($id_str)";

     #   my $fs = $self->generic_fetch($constraint, $mapper, $slice);

 #       $fs = _remap($fs, $mapper, $slice);

  #      push @features, @$fs;

   #   } else {
        # do multiple split queries using start / end constraints

    #    my $max_len = $self->_max_feature_length() ||
      #    $mcc->fetch_max_length_by_CoordSystem_feature_type($feat_cs,
    #                                                         $tab_name);
     #   my $len = @coords;
       # for(my $i = 0; $i < $len; $i++) {
       #   my $constraint = $orig_constraint;
       #   $constraint .= " AND " if($constraint);
       #   $constraint .=
       #       "${tab_syn}.seq_region_id = "     . $ids[$i] . " AND " .
       #       "${tab_syn}.seq_region_start <= " . $coords[$i]->end() . " AND ".
       #       "${tab_syn}.seq_region_end >= "   . $coords[$i]->start();

       #   if($max_len) {
       #     my $min_start = $coords[$i]->start() - $max_len;
       #     $constraint .=
       #       " AND ${tab_syn}.seq_region_start >= $min_start";
       #   }

        #  my $fs = $self->generic_fetch($constraint,$mapper,$slice);

 #         $fs = _remap($fs, $mapper, $slice);

#          push @features, @$fs;
    #    }
    #  }
  ###  }
  } #COORD system loop

  return \@features;
}




#
# Given a list of features checks if they are in the correct coord system
# by looking at the first features slice.  If they are not then they are
# converted and placed on the slice.
#
sub _remap {
  my ($features, $mapper, $slice) = @_;

  #check if any remapping is actually needed
  if(@$features && (!$features->[0]->isa('Bio::EnsEMBL::Feature') ||
                    $features->[0]->slice == $slice)) {
    return $features;
  }

  #remapping has not been done, we have to do our own conversion from
  #to slice coords

  my @out;

  my $slice_start = $slice->start();
  my $slice_end   = $slice->end();
  my $slice_strand = $slice->strand();
  my $slice_cs    = $slice->coord_system();

  my ($seq_region, $start, $end, $strand);

  #my $slice_seq_region_id = $slice->get_seq_region_id();
  my $slice_seq_region = $slice->seq_region_name();

  foreach my $f (@$features) {
    #since feats were obtained in contig coords, attached seq is a contig
    my $fslice = $f->slice();
    if(!$fslice) {
      throw("Feature does not have attached slice.\n");
    }
    my $fseq_region = $fslice->seq_region_name();
    my $fseq_region_id = $fslice->get_seq_region_id();
    my $fcs = $fslice->coord_system();

    if(!$slice_cs->equals($fcs)) {
      #slice of feature in different coord system, mapping required

      ($seq_region, $start, $end, $strand) =
        $mapper->fastmap($fseq_region_id,$f->start(),$f->end(),$f->strand(),$fcs);

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



=head2 fetch_all_by_external_name

  Arg [1]    : String $external_name
               An external identifier of the feature to be obtained
  Arg [2]    : (optional) String $external_db_name
               The name of the external database from which the
               identifier originates.
  Example    : my @features =
                  @{ $adaptor->fetch_all_by_external_name( 'NP_065811.1') };
  Description: Retrieves all features which are associated with
               an external identifier such as a GO term, Swissprot
               identifer, etc.  Usually there will only be a single
               feature returned in the list reference, but not
               always.  Features are returned in their native
               coordinate system, i.e. the coordinate system in which
               they are stored in the database.  If they are required
               in another coordinate system the Feature::transfer or
               Feature::transform method can be used to convert them.
               If no features with the external identifier are found,
               a reference to an empty list is returned.
  Returntype : arrayref of Bio::EnsEMBL::Feature objects
  Exceptions : none
  Caller     : general
  Status     : at risk

=cut

sub fetch_all_by_external_name {
  my ( $self, $external_name, $external_db_name ) = @_;

  my $entryAdaptor = $self->db->get_DBEntryAdaptor();

  my @ids;

  if($self->isa('Bio::EnsEMBL::Funcgen::DBSQL::RegulatoryFeatureAdaptor')){
	@ids = $entryAdaptor->list_regulatory_feature_ids_by_extid($external_name, $external_db_name);
  }else{
	warn "Does not yet accomodate non regulatory feature external names";
	return;
	@ids = $entryAdaptor->list_external_feature_ids_by_extid($external_name, $external_db_name);
  }

  return $self->fetch_all_by_dbID_list( \@ids );
}

=head2 fetch_by_display_label

  Arg [1]    : String $label - display label of feature to fetch
  Example    : my $feat = $adaptor->fetch_by_display_label("BRCA2");
  Description: Returns the feature which has the given display label or
               undef if there is none. If there are more than 1, only the first
               is reported.
  Returntype : Bio::EnsEMBL::Funcgen::Feature
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut

sub fetch_by_display_label {
  my $self = shift;
  my $label = shift;

  my @tables = $self->_tables();
  #my $syn = $tables[0]->[1];

  my $constraint = "x.display_label = '$label'";# AND ${syn}.is_current = 1";
  my ($feature) = @{ $self->generic_fetch($constraint) };

  return $feature;
}





1;


