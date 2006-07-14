#
# Ensembl module for Bio::EnsEMBL::DBSQL::Funcgen::OligoFeatureAdaptor
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::DBSQL::Funcgen::OligoFeatureAdaptor - A database adaptor for fetching and
storing OligoFeature objects.

=head1 SYNOPSIS

my $ofa = $db->get_OligoFeatureAdaptor();

my $features = $ofa->fetch_all_by_Probe($probe);
$features = $ofa->fetch_all_by_Slice_arrayname($slice, 'Array-1', 'Array-2');

=head1 DESCRIPTION

The OligoFeatureAdaptor is a database adaptor for storing and retrieving
OligoFeature objects.

=head1 AUTHOR

This module was created by Ian Sealy, but is almost entirely based on the
OligoFeatureAdaptor module written by Arne Stabenau.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::DBSQL::OligoFeatureAdaptor;

use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Funcgen::OligoFeature;
use Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor);


=head2 fetch_all_by_Probe

  Arg [1]    : Bio::EnsEMBL::Funcgen::OligoProbe
  Example    : my $features = $ofa->fetch_all_by_Probe($probe);
  Description: Fetchs all features that a given probe creates.
  Returntype : Listref of Bio::EnsEMBL::OligoFeature objects
  Exceptions : Throws if argument is not a stored OligoProbe object
  Caller     : OligoProbe->get_all_OligoFeatures()
  Status     : Medium Risk

=cut

sub fetch_all_by_Probe {
	my $self  = shift;
	my $probe = shift;
	
	if ( !ref($probe) && !$probe->isa('Bio::EnsEMBL::Funcgen::OligoProbe') ) {
		throw('fetch_all_by_Probe requires a Bio::EnsEMBL::Funcgen::OligoProbe object');
	}
	if ( !defined $probe->dbID() ) {
		throw('fetch_all_by_Probe requires a stored Bio::EnsEMBL::Funcgen::OligoProbe object');
	}
	
	return $self->generic_fetch( 'of.oligo_probe_id = ' . $probe->dbID() );
}

=head2 fetch_all_by_probeset

  Arg [1]    : string - probeset
  Example    : my $features = $ofa->fetch_all_by_probeset('Set-1');
  Description: Fetchs all features that a given probeset creates.
  Returntype : Listref of Bio::EnsEMBL::OligoFeature objects
  Exceptions : Throws if no probeset argument
  Caller     : General
  Status     : Medium Risk

=cut

sub fetch_all_by_probeset {
	my $self     = shift;
	my $probeset = shift;
	
	throw("Not implmeneted\n");

	if (!$probeset) {
		throw('fetch_all_by_probeset requires a probeset argument');
	}
	
	

	return $self->generic_fetch( "op.probeset = '$probeset'" );
}


#Need to add:
#fetch_all_by_Slice_Experiment
#fetch_all_by_Slice_experimentname ? name not unique enough?


=head2 fetch_all_by_Slice_arrayname

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2...] : List of strings - array name(s)
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $features = $ofa->fetch_by_Slice_arrayname($slice, '');
  Description: Retrieves a list of features on a given slice that are created
               by probes from the specified arrays.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::OligoFeature objects
  Exceptions : Throws if no array name is provided
  Caller     : Slice->get_all_OligoFeatures()
  Status     : Medium Risk

=cut

sub fetch_all_by_Slice_arrayname {
	my ($self, $slice, @arraynames) = @_;
	
	throw('Need array name as parameter') if !@arraynames;
	
	my $constraint;
	if (scalar @arraynames == 1) {
		#Will this work
		#will this pick up the array_chip_id link from array_chip to oligo_probe?
		$constraint = qq( a.name = '$arraynames[0]' AND a.array_id = ac.array_id );
		#$constraint = qq( a.name = '$arraynames[0]' );
	} else {
		throw("Not implemented for multple arrays");
		$constraint = join q(','), @arraynames;
		$constraint = qq( a.name IN ('$constraint') );
	}
	
	return $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint);
}

=head2 fetch_all_by_Slice_type

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : string - type of array (e.g. AFFY or OLIGO)
  Arg [3]    : (optional) string - logic name
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $features = $ofa->fetch_by_Slice_type($slice, 'OLIGO');
  Description: Retrieves a list of features on a given slice that are created
               by probes from the specified type of array.
  Returntype : Listref of Bio::EnsEMBL::OligoFeature objects
  Exceptions : Throws if no array type is provided
  Caller     : General
  Status     : Medium Risk

=cut

sub fetch_all_by_Slice_type {
	my ($self, $slice, $type, $logic_name) = @_;

	throw("Not implemented yet\n");
	
	throw('Need type as parameter') if !$type;
	
	my $constraint = qq( a.type = '$type' );
	
	return $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint, $logic_name);
}
 
=head2 _tables

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns the names and aliases of the tables to use for queries.
  Returntype : List of listrefs of strings
  Exceptions : None
  Caller     : Internal
  Status     : Medium Risk

=cut

sub _tables {
	my $self = shift;
	
	return (
			[ 'oligo_feature', 'of' ], 
			[ 'oligo_probe',   'op' ], 
			[ 'array_chip',   'ac' ],
			[ 'array',   'a' ]
		   );
}

=head2 _columns

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns a list of columns to use for queries.
  Returntype : List of strings
  Exceptions : None
  Caller     : Internal
  Status     : Medium Risk

=cut

sub _columns {
	my $self = shift;
	
	return qw(
		of.oligo_feature_id  of.seq_region_id
		of.seq_region_start  of.seq_region_end
		of.seq_region_strand of.mismatches
		of.oligo_probe_id    of.analysis_id
		a.name              op.probeset
		op.name
	);
}

=head2 _default_where_clause

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns an additional table joining constraint to use for
			   queries.
  Returntype : List of strings
  Exceptions : None
  Caller     : Internal
  Status     : Medium Risk

=cut
sub _default_where_clause {
	my $self = shift;
	
	return 'of.oligo_probe_id = op.oligo_probe_id AND op.array_chip_id = ac.array_chip_id';
}

=head2 _final_clause

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns an ORDER BY clause. Sorting by oligo_feature_id would be
			   enough to eliminate duplicates, but sorting by location might
			   make fetching features on a slice faster.
  Returntype : String
  Exceptions : None
  Caller     : generic_fetch
  Status     : Medium Risk

=cut


sub _final_clause {
	return ' ORDER BY of.seq_region_id, of.seq_region_start, of.oligo_feature_id';
}


=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates OligoFeature objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::OligoFeature objects
  Exceptions : None
  Caller     : Internal
  Status     : Medium Risk

=cut

sub _objs_from_sth {
	my ($self, $sth, $mapper, $dest_slice) = @_;
	
	# This code is ugly because caching is used to improve speed

	my $sa = $self->db->get_SliceAdaptor();
	my $aa = $self->db->get_AnalysisAdaptor();

	my @features;
	
	my (%analysis_hash, %slice_hash, %sr_name_hash, %sr_cs_hash);

	my (
		$oligo_feature_id,  $seq_region_id,
		$seq_region_start,  $seq_region_end,
		$seq_region_strand, $mismatches,
		$oligo_probe_id,    $analysis_id,
		$array_name,        $probeset,
		$oligo_probe_name,
	);
	$sth->bind_columns(
		\$oligo_feature_id,  \$seq_region_id,
		\$seq_region_start,  \$seq_region_end,
		\$seq_region_strand, \$mismatches,
		\$oligo_probe_id,    \$analysis_id,
		\$array_name,        \$probeset,
		\$oligo_probe_name,
	);

	my $asm_cs;
	my $cmp_cs;
	my $asm_cs_name;
	my $asm_cs_vers;
	my $cmp_cs_name;
	my $cmp_cs_vers;
	if ($mapper) {
		$asm_cs      = $mapper->assembled_CoordSystem();
		$cmp_cs      = $mapper->component_CoordSystem();
		$asm_cs_name = $asm_cs->name();
		$asm_cs_vers = $asm_cs->version();
		$cmp_cs_name = $cmp_cs->name();
		$cmp_cs_vers = $cmp_cs->version();
	}

	my $dest_slice_start;
	my $dest_slice_end;
	my $dest_slice_strand;
	my $dest_slice_length;
	my $dest_slice_sr_name;
	if ($dest_slice) {
		$dest_slice_start   = $dest_slice->start();
		$dest_slice_end     = $dest_slice->end();
		$dest_slice_strand  = $dest_slice->strand();
		$dest_slice_length  = $dest_slice->length();
		$dest_slice_sr_name = $dest_slice->seq_region_name();
	}

	my $last_feature_id = -1;
	FEATURE: while ( $sth->fetch() ) {

		# This assumes that features come out sorted by ID
		next if ($last_feature_id == $oligo_feature_id);
		$last_feature_id = $oligo_feature_id;

		# Get the analysis object
		my $analysis = $analysis_hash{$analysis_id} ||= $aa->fetch_by_dbID($analysis_id);

		# Get the slice object
		my $slice = $slice_hash{'ID:'.$seq_region_id};

		if (!$slice) {
			$slice                            = $sa->fetch_by_seq_region_id($seq_region_id);
			$slice_hash{'ID:'.$seq_region_id} = $slice;
			$sr_name_hash{$seq_region_id}     = $slice->seq_region_name();
			$sr_cs_hash{$seq_region_id}       = $slice->coord_system();
		}

		my $sr_name = $sr_name_hash{$seq_region_id};
		my $sr_cs   = $sr_cs_hash{$seq_region_id};

		# Remap the feature coordinates to another coord system if a mapper was provided
		if ($mapper) {
			($sr_name, $seq_region_start, $seq_region_end, $seq_region_strand)
				= $mapper->fastmap($sr_name, $seq_region_start, $seq_region_end, $seq_region_strand, $sr_cs);

			# Skip features that map to gaps or coord system boundaries
			next FEATURE if !defined $sr_name;

			# Get a slice in the coord system we just mapped to
			if ( $asm_cs == $sr_cs || ( $cmp_cs != $sr_cs && $asm_cs->equals($sr_cs) ) ) {
				$slice = $slice_hash{"NAME:$sr_name:$cmp_cs_name:$cmp_cs_vers"}
					||= $sa->fetch_by_region($cmp_cs_name, $sr_name, undef, undef, undef, $cmp_cs_vers);
			} else {
				$slice = $slice_hash{"NAME:$sr_name:$asm_cs_name:$asm_cs_vers"}
					||= $sa->fetch_by_region($asm_cs_name, $sr_name, undef, undef, undef, $asm_cs_vers);
			}
		}

		# If a destination slice was provided convert the coords
		# If the destination slice starts at 1 and is forward strand, nothing needs doing
		if ($dest_slice) {
			unless ($dest_slice_start == 1 && $dest_slice_strand == 1) {
				if ($dest_slice_strand == 1) {
					$seq_region_start = $seq_region_start - $dest_slice_start + 1;
					$seq_region_end   = $seq_region_end   - $dest_slice_start + 1;
				} else {
					my $tmp_seq_region_start = $seq_region_start;
					$seq_region_start        = $dest_slice_end - $seq_region_end       + 1;
					$seq_region_end          = $dest_slice_end - $tmp_seq_region_start + 1;
					$seq_region_strand      *= -1;
				}
			}

			# Throw away features off the end of the requested slice
			next FEATURE if $seq_region_end < 1 || $seq_region_start > $dest_slice_length
				|| ( $dest_slice_sr_name ne $sr_name );

			$slice = $dest_slice;
		}

		push @features, $self->_new_fast( {
			'start'         => $seq_region_start,
			'end'           => $seq_region_end,
			'strand'        => $seq_region_strand,
			'slice'         => $slice,
			'analysis'      => $analysis,
			'adaptor'       => $self,
			'dbID'          => $oligo_feature_id,
			'mismatchcount' => $mismatches,
			'_probe_id'     => $oligo_probe_id,
			'probeset'      => $probeset,
			'_probe_name'   => $oligo_probe_name
		} );
	}

	return \@features;
}

=head2 _new_fast

  Args       : Hashref to be passed to OligoFeature->new_fast()
  Example    : None
  Description: Construct an OligoFeature object using quick and dirty new_fast.
  Returntype : Bio::EnsEMBL::Funcgen::OligoFeature
  Exceptions : None
  Caller     : _objs_from_sth
  Status     : Medium Risk

=cut

sub _new_fast {
	my $self = shift;
	
	my $hash_ref = shift;
	return Bio::EnsEMBL::Funcgen::OligoFeature->new_fast($hash_ref);
}

=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::OligoFeature objects
  Example    : $ofa->store(@features);
  Description: Stores given OligoFeature objects in the database. Should only
               be called once per feature because no checks are made for
			   duplicates. Sets dbID and adaptor on the objects that it stores.
  Returntype : None
  Exceptions : Throws if a list of OligoFeature objects is not provided or if
               an analysis is not attached to any of the objects
  Caller     : General
  Status     : Medium Risk

=cut

sub store{
	my ($self, @ofs) = @_;

	if (scalar(@ofs) == 0) {
		throw('Must call store with a list of OligoFeature objects');
	}

	my $sth = $self->prepare("
		INSERT INTO oligo_feature (
			seq_region_id,  seq_region_start,
			seq_region_end, seq_region_strand,
			oligo_probe_id,  analysis_id,
			mismatches
		) VALUES (?, ?, ?, ?, ?, ?, ?)
	");

	my $db = $self->db();
	my $analysis_adaptor = $db->get_AnalysisAdaptor();

	FEATURE: foreach my $of (@ofs) {

		if( !ref $of || !$of->isa('Bio::EnsEMBL::Funcgen::OligoFeature') ) {
			throw('Feature must be an OligoFeature object');
		}

		if ( $of->is_stored($db) ) {
			warning('OligoFeature [' . $of->dbID() . '] is already stored in the database');
			next FEATURE;
		}

		if ( !defined $of->analysis() ) {
			throw('An analysis must be attached to the OligoFeature objects to be stored.');
		}

		# Store the analysis if it has not been stored yet
		if ( !$of->analysis->is_stored($db) ) {
			$analysis_adaptor->store( $of->analysis() );
		}

		my $original = $of;
		my $seq_region_id;
		($of, $seq_region_id) = $self->_pre_store($of);

		$sth->bind_param(1, $seq_region_id,        SQL_INTEGER);
		$sth->bind_param(2, $of->start(),          SQL_INTEGER);
		$sth->bind_param(3, $of->end(),            SQL_INTEGER);
		$sth->bind_param(4, $of->strand(),         SQL_TINYINT);
		$sth->bind_param(5, $of->probe->dbID(),    SQL_INTEGER);
		$sth->bind_param(6, $of->analysis->dbID(), SQL_INTEGER);
		$sth->bind_param(7, $of->mismatchcount(),  SQL_TINYINT);

		$sth->execute();

		$original->dbID( $sth->{'mysql_insertid'} );
		$original->adaptor($self);
	}
}

=head2 list_dbIDs

  Args       : None
  Example    : my @feature_ids = @{$ofa->list_dbIDs()};
  Description: Gets an array of internal IDs for all OligoFeature objects in
               the current database.
  Returntype : List of ints
  Exceptions : None
  Caller     : ?
  Status     : Medium Risk

=cut

sub list_dbIDs {
	my $self = shift;
	
	return $self->_list_dbIDs('oligo_feature');
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

  #my $slice_adaptor = $db->get_SliceAdaptor();
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


  #Need to add to Funcgen coord_system here
  #check if name and version are present and reset coord_system_id to that one, else get last ID and create a new one
  #coord_system_ids will not match those in core DBs, so we need ot be mindful about this.
  #can't use is_stored as this simply checks the dbID
  #seq_region_ids may change between shemas with the same assembly version
  #Need to store schema_version somewhere to maintain the seq_region_id mapping
  #extend the coord_system table to have schema version, link to feature tables via coord_system_id
  #how are we going to retrieve, we have species and schema version, so will have to do some jiggery pokery
  #There is a possibility that the same schema may be used for two different gene builds
  #thus having the same name version schema triplets, but different seq_region_ids
  #can't really code around this...unless we stored the schema and the build instead of just the schema
  #this would make it easier to generate the dnadb as we could simply concat $species."_core_".$schema_build
  #can not get this from the meta table, so we'll have to fudge it from the db_name
  
  #Do we need to check the the dnadb and the slice db match?
  #Do we have to have specified a dnadb at this point?  No.
  #But need to put checks in place for dnadb methods i.e. seq/slice retrieval

  #Build object cache here to prevent multiple(possible 1000's of calls to the coord tables;
  $self->{'_ofa_coord_cache'}  ||= {};
  #a bit hack this, can we not make dnadb mandatory and set in new?!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #would still have to do tests to check we're not passing a slice from a db which isn't the dnadb :(

  #This will get called for each feature!!
  #No guarantee that each feature will be built from the same db (could set once in store otherwise)
  #my $schema_build = ${$slice->adaptor->db->db_handle->selectrow_array("SELECT meta_value value from meta where mate_key = \"data.version\"")}[0];  #do we need to check whether there is only one?
  my $schema_build = $slice->adaptor->db->db_name();
  $schema_build =~ s/[a-zA-Z_]//;
  print "schema_build id $schema_build\n";

  if(! exists $self->{'_ofa_coord_cache'}{$schema_build.":".$cs->version().":".$cs->name()}){
	  #This should now only be called for the number of unique name version schema triplets
	  #In reality this will only be unique names for the given version we are creating the feature on
	  #e.g. all the chromosomes
	  my $csa = $self->db->get_CoordSystemAdaptor();
	
	  #store coordsystem with schema_version
	  #This will enable automatic dnadb DBAdaptor generation
	  #otherwise would have to list all DBs in registry and select the one which matched the schema version
	  #This will also enable validation if a non-standard dnadb is passed for retrieval
	  #can check meta table for schema.version (data.version? genebuild.name?)

	  #Generate Funcgen::CoordSystem
	  $cs = Bio::EnsEMBL::Funcgen::CoordSystem->new(
													-NAME    => $cs->name(),
                                                    -VERSION => $cs->version(),
                                                    -RANK    => $cs->rank(),
                                                    #-ADAPTOR => $adaptor,
                                                    -DEFAULT => $cs->default(),
                                                    -SEQUENCE_LEVEL => $cs->is_sequence_level(),
													-SCHEMA_BUILD => $schema_build
												   );

	  $csa->store($cs);

  }


  #retrieve correpsonding Funcgen coord_system
  $cs = $csa->fetch_by_name_version_schema_build($cs->name(), $cs->version(), $schema_build);

  #Funcgen::CoordSystemAdaptor
  #fetch_by_name_version_schema_build
  #remove all other fetch's and innapropriate methods?
  
  #Funcgen::CoordSystem
  #add methods to retrieve data from dnadb/core coord_system tables
  

  my ($tab) = $self->_tables();
  my $tabname = $tab->[0];


  #Need to do this for Funcgen DB
  my $mcc = $db->get_MetaCoordContainer();
  $mcc->add_feature_type($cs, $tabname, $feature->length);

 # my $seq_region_id = $slice_adaptor->get_seq_region_id($slice);
  my $seq_region_id = $slice->get_seq_region_id();

  #would never get called as we're not validating against a different core DB
  #if(!$seq_region_id) {
  #  throw('Feature is associated with seq_region which is not in this dnadb.');
  #}

  return ($feature, $seq_region_id);
}



1;

