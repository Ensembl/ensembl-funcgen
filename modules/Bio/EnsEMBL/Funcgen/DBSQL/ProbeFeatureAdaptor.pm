#
# Ensembl module for Bio::EnsEMBL::DBSQL::Funcgen::ProbeFeatureAdaptor
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::DBSQL::Funcgen::ProbeFeatureAdaptor - A database adaptor for fetching and
storing ProbeFeature objects.

=head1 SYNOPSIS

my $ofa = $db->get_ProbeFeatureAdaptor();

my $features = $ofa->fetch_all_by_Probe($probe);
$features = $ofa->fetch_all_by_Slice_arrayname($slice, 'Array-1', 'Array-2');

=head1 DESCRIPTION

The ProbeFeatureAdaptor is a database adaptor for storing and retrieving
ProbeFeature objects.

=head1 AUTHOR

This module was created by Nathan Johnson.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::DBSQL::ProbeFeatureAdaptor;

use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Funcgen::ProbeFeature;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;
use vars qw(@ISA);
use strict;
use warnings;

@ISA = qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor);


=head2 fetch_all_by_Probe

  Arg [1]    : Bio::EnsEMBL::Funcgen::Probe
  Example    : my $features = $ofa->fetch_all_by_Probe($probe);
  Description: Fetchs all features that a given probe creates.
  Returntype : Listref of Bio::EnsEMBL::PasteFeature objects
  Exceptions : Throws if argument is not a stored Probe object
  Caller     : Probe->get_all_ProbeFeatures()
  Status     : At Risk

=cut

sub fetch_all_by_Probe {
  my $self  = shift;
  my $probe = shift;
  
  if ( !ref($probe) && !$probe->isa('Bio::EnsEMBL::Funcgen::Probe') ) {
    throw('fetch_all_by_Probe requires a Bio::EnsEMBL::Funcgen::Probe object');
  }

  if ( !defined $probe->dbID() ) {
    throw('fetch_all_by_Probe requires a stored Bio::EnsEMBL::Funcgen::Probe object');
  }
	
  return $self->generic_fetch( 'pf.probe_id = ' . $probe->dbID() );
}

=head2 fetch_all_by_Probe_id

  Arg [1]    : int - Probe dbID
  Example    : my @features = @{$ofa->fetch_all_by_Probe_id($pid)};
  Description: Fetchs all features that a given probe creates.
  Returntype : Listref of Bio::EnsEMBL::PasteFeature objects
  Exceptions : Throws if argument not defined
  Caller     : Probe->get_all_ProbeFeatures()
  Status     : At Risk

=cut

sub fetch_all_by_Probe_id {
  my $self  = shift;
  my $pid = shift;
  
  if ( ! defined $pid ) {
    throw('Need to specify a probe _id');
  }
	
  return $self->generic_fetch( 'pf.probe_id = ' . $pid );
}



=head2 fetch_all_by_probeset

  Arg [1]    : string - probeset
  Example    : my $features = $ofa->fetch_all_by_probeset('Set-1');
  Description: Fetchs all features that a given probeset creates.
  Returntype : Listref of Bio::EnsEMBL::ProbeFeature objects
  Exceptions : Throws if no probeset argument
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_probeset {
	my $self     = shift;
	my $probeset = shift;
	
	throw("Not implmeneted\n");

	if (!$probeset) {
		throw('fetch_all_by_probeset requires a probeset argument');
	}
	
	

	return $self->generic_fetch( "p.probeset = '$probeset'" );
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
  Returntype : Listref of Bio::EnsEMBL::Funcgen::ProbeFeature objects
  Exceptions : Throws if no array name is provided
  Caller     : Slice->get_all_ProbesFeatures()
  Status     : At Risk

=cut

sub fetch_all_by_Slice_arrayname {
	my ($self, $slice, @arraynames) = @_;

	throw("This should return data from all experiments, but will break if arrays are mapped to different coord_systems");
	
	throw('Need array name as parameter') if !@arraynames;
	
	my $constraint;
	if (scalar @arraynames == 1) {
		#Will this work
		#will this pick up the array_chip_id link from array_chip to probe?
		$constraint = qq( a.name = '$arraynames[0]' AND a.array_id = ac.array_id );
		#$constraint = qq( a.name = '$arraynames[0]' );
	} else {
		throw("Not implemented for multple arrays");
		$constraint = join q(','), @arraynames;
		$constraint = qq( a.name IN ('$constraint') );
	}
	
	return $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint);
}


#should this take >1 EC? What if we can't fit a all mappings onto one chip
#Would possibly miss some from the slice

=head2 fetch_all_by_Slice_ExperimentalChips

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2...] : listref of Bio::EnsEMBL::Funcgen::ExperimentalChip objects
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $features = $ofa->fetch_by_Slice_arrayname($slice, $exp);
  Description: Retrieves a list of features on a given slice that are created
               by probes from the given ExperimentalChip.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::ProbeFeature objects
  Exceptions : Throws if no array name is provided
  Caller     : 
  Status     : At Risk

=cut

sub fetch_all_by_Slice_ExperimentalChips {
  my ($self, $slice, $exp_chips) = @_;

  my (%nr);


  foreach my $ec(@$exp_chips){
    
    throw("Need pass listref of valid Bio::EnsEMBL::Funcgen::ExperimentalChip objects") 
      if ! $ec->isa("Bio::EnsEMBL::Funcgen::ExperimentalChip");
    
    $nr{$ec->array_chip_id()} = 1;
  }
   
  my $constraint = " p.array_chip_id IN (".join(", ", keys %nr).") AND p.probe_id = pf.probe_id ";

    
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
  Returntype : Listref of Bio::EnsEMBL::ProbeFeature objects
  Exceptions : Throws if no array type is provided
  Caller     : General
  Status     : At Risk

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
  Status     : At Risk

=cut

sub _tables {
	my $self = shift;
	
	return (
			[ 'probe_feature', 'pf' ], 
			[ 'probe',   'p' ], 
			#[ 'array_chip',   'ac' ],#these are required for array based queries not implemented yet
			#[ 'array',   'a' ]
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
  Status     : At Risk

=cut

sub _columns {
	my $self = shift;


	#do we need array_name?
	
	return qw(
			  pf.probe_feature_id  pf.seq_region_id
			  pf.seq_region_start  pf.seq_region_end
			  pf.seq_region_strand pf.probe_id    
			  pf.analysis_id	   pf.mismatches
			  pf.cigar_line        p.name
			 );

	#removed probeset and array name
	
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
  Status     : At Risk

=cut
sub _default_where_clause {
	my $self = shift;
	
	return 'pf.probe_id = p.probe_id';# AND p.array_chip_id = ac.array_chip_id';
}

=head2 _final_clause

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns an ORDER BY clause. Sorting by probe_feature_id would be
			   enough to eliminate duplicates, but sorting by location might
			   make fetching features on a slice faster.
  Returntype : String
  Exceptions : None
  Caller     : generic_fetch
  Status     : At Risk

=cut


sub _final_clause {
	return ' ORDER BY pf.seq_region_id, pf.seq_region_start, pf.probe_feature_id';
}


=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates ProbeFeature objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::ProbeFeature objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
	my ($self, $sth, $mapper, $dest_slice) = @_;

	#For EFG this has to use a dest_slice from core/dnaDB whether specified or not.
	#So if it not defined then we need to generate one derived from the species_name and schema_build of the feature we're retrieving.


	
	# This code is ugly because caching is used to improve speed
	my ($sa, $seq_region_id);
	$sa = $dest_slice->adaptor->db->get_SliceAdaptor() if($dest_slice);#don't really need this if we're using DNADBSliceAdaptor?

	#Some of this in now probably overkill as we'll always be using the DNADB as the slice DB
	#Hence it should always be on the same coord system

	my $aa = $self->db->get_AnalysisAdaptor();
	my @features;
	my (%analysis_hash, %slice_hash, %sr_name_hash, %sr_cs_hash);

	my (
	    $probe_feature_id,  $efg_seq_region_id,
	    $seq_region_start,  $seq_region_end,
	    $seq_region_strand, $mismatches,
		$probe_id,    	    $analysis_id,
		$probe_name,	    $cigar_line,
	);
	$sth->bind_columns(
					   \$probe_feature_id,  \$efg_seq_region_id,
					   \$seq_region_start,  \$seq_region_end,
					   \$seq_region_strand, \$probe_id,
					   \$analysis_id,       \$mismatches,
					   \$cigar_line,        \$probe_name
	);

	my ($asm_cs, $cmp_cs, $asm_cs_name, $asm_cs_vers ,$cmp_cs_name, $cmp_cs_vers);
	if ($mapper) {
		$asm_cs      = $mapper->assembled_CoordSystem();
		$cmp_cs      = $mapper->component_CoordSystem();
		$asm_cs_name = $asm_cs->name();
		$asm_cs_vers = $asm_cs->version();
		$cmp_cs_name = $cmp_cs->name();
		$cmp_cs_vers = $cmp_cs->version();
	}

	my ($dest_slice_start, $dest_slice_end, $dest_slice_strand);
	my ($dest_slice_length, $dest_slice_sr_name);
	if ($dest_slice) {
		$dest_slice_start   = $dest_slice->start();
		$dest_slice_end     = $dest_slice->end();
		$dest_slice_strand  = $dest_slice->strand();
		$dest_slice_length  = $dest_slice->length();
		$dest_slice_sr_name = $dest_slice->seq_region_name();
	}

	#This has already been done by
	#build seq_region_cache based on slice
	#$self->build_seq_region_cache_by_Slice($slice);

	FEATURE: while ( $sth->fetch() ) {
		  #Need to build a slice adaptor cache here?
		  #Would only ever want to do this if we enable mapping between assemblies??
		  #Or if we supported the mapping between cs systems for a given schema_build, which would have to be handled by the core api
		  
		#get core seq_region_id
		$seq_region_id = $self->get_core_seq_region_id($efg_seq_region_id);
		
		if(! $seq_region_id){
		  warn "Cannot get slice for eFG seq_region_id $efg_seq_region_id\n".
		  "The region you are using is not present in the cuirrent dna DB";
		  next;
		}

		
		#warn "Need to implement slice adaptor hash, based on seq_region id??";#


		#we need to be mindful of dynamic assembly mapping
		#or would this be handled before here?
		#will different cs_id be handled before here also, so we would never see different cs_ids?

		#if($old_cs_id && ($old_cs_id != $cs_id)){
		#  throw("More than one coord_system for feature query, need to implement SliceAdaptor hash?");
		#}
		
		#$old_cs_id = $cs_id;


		#This should by default be the slice adaptor of the dnadb we're concerned with
		#what about assembly mapping where we have feature from multiple dnadbs returned in the same query


		#this needs to be reset for each seq_region_id
		$sa ||= $self->db->get_SliceAdaptor();#$cs_id);


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

		#need to check once more here as it may not be in the DB, 
		#i.e. a supercontig(non-versioned) may have been deleted between releases


		my $sr_name = $sr_name_hash{$seq_region_id};
		my $sr_cs   = $sr_cs_hash{$seq_region_id};

		# Remap the feature coordinates to another coord system if a mapper was provided
		if ($mapper) {

			throw("Not yet implmented mapper, check equals are Funcgen calls too!");

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


		warn "Got ".scalar(@features)." features" if (! (scalar(@features)%1000));

		push @features, $self->_new_fast( {
											 'start'         => $seq_region_start,
											 'end'           => $seq_region_end,
											 'strand'        => $seq_region_strand,
											 'slice'         => $slice,
											 'analysis'      => $analysis,#we should lazy load this
											 'adaptor'       => $self,
											 'dbID'          => $probe_feature_id,
											 'mismatchcount' => $mismatches,
											 'cigar_line'    => $cigar_line,
											 'probe_id'     => $probe_id,
											 #'probeset'      => $probeset,#???do we need this?
											 '_probe_name'   => $probe_name
											} );


	
	  }

	return \@features;
}

=head2 _new_fast

  Args       : Hashref to be passed to ProbeFeature->new_fast()
  Example    : None
  Description: Construct an ProbeFeature object using quick and dirty new_fast.
  Returntype : Bio::EnsEMBL::Funcgen::ProbeFeature
  Exceptions : None
  Caller     : _objs_from_sth
  Status     : Medium Risk

=cut

sub _new_fast {
	my $self = shift;
	
	my $hash_ref = shift;
	return Bio::EnsEMBL::Funcgen::ProbeFeature->new_fast($hash_ref);
}

=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::ProbeFeature objects
  Example    : $ofa->store(@features);
  Description: Stores given ProbeFeature objects in the database. Should only
               be called once per feature because no checks are made for
			   duplicates. Sets dbID and adaptor on the objects that it stores.
  Returntype : None
  Exceptions : Throws if a list of ProbeFeature objects is not provided or if
               an analysis is not attached to any of the objects
  Caller     : General
  Status     : At Risk

=cut

sub store{
	my ($self, @ofs) = @_;

	if (scalar(@ofs) == 0) {
		throw('Must call store with a list of ProbeFeature objects');
	}

	my $sth = $self->prepare("
		INSERT INTO probe_feature (
			seq_region_id,  seq_region_start,
			seq_region_end, seq_region_strand,
          	probe_id,  analysis_id,
			mismatches, cigar_line
		) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
	");

	my $db = $self->db();
	my $analysis_adaptor = $db->get_AnalysisAdaptor();

	FEATURE: foreach my $of (@ofs) {

		if( !ref $of || !$of->isa('Bio::EnsEMBL::Funcgen::ProbeFeature') ) {
			throw('Feature must be an ProbeFeature object');
		}

		if ( $of->is_stored($db) ) {
			warning('ProbeFeature [' . $of->dbID() . '] is already stored in the database');
			next FEATURE;
		}

		if ( !defined $of->analysis() ) {
			throw('An analysis must be attached to the ProbeFeature objects to be stored.');
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
		$sth->bind_param(5, $of->probe_id(),      SQL_INTEGER);
		$sth->bind_param(6, $of->analysis->dbID(), SQL_INTEGER);
		$sth->bind_param(7, $of->mismatchcount(),  SQL_TINYINT);
		$sth->bind_param(8, $of->cigar_line(),     SQL_VARCHAR);

		$sth->execute();

		$original->dbID( $sth->{'mysql_insertid'} );
		$original->adaptor($self);
	}

	return \@ofs
}

=head2 list_dbIDs

  Args       : None
  Example    : my @feature_ids = @{$ofa->list_dbIDs()};
  Description: Gets an array of internal IDs for all ProbeFeature objects in
               the current database.
  Returntype : List of ints
  Exceptions : None
  Caller     : ?
  Status     : Medium Risk

=cut

sub list_dbIDs {
	my $self = shift;
	
	return $self->_list_dbIDs('probe_feature');
}


=head2 reassign_features_to_probe

  Arg[0]     : ARRAYREF - feature dbIDs to reassign
  Arg[1]     : int - probe dbID to reassign to
  Example    : $ofa->reassign_feature_to_probe(\@fids, $pid);
  Description: Update features to link to given probe dbID
  Returntype : None
  Exceptions : Throws is args not met
  Caller     : Importer
  Status     : At Risk

=cut

sub reassign_feature_to_probe{
	my ($self, $fids_ref, $pid) = @_;
	
	if(! @$fids_ref || ! $pid){
	  throw('Need to pass a ref to an array of feature ids and a probe id to reassign to');
	}
	
	my $cmd = 'UPDATE probe_feature SET probe_id='.$pid.' WHERE probe_feature_id IN ('.join(',', @$fids_ref).')';
	$self->db->dbc->do($cmd);

	#This will fail anyway?
	#if($?){
	#  throw("SQL Command failed:\t$sql\n$@");
	#}

	return;
}

=head2 delete_features

  Arg[0]     : ARRAYREF - feature dbIDs to reassign
  Example    : $pfa->delete_feature(\@fids);
  Description: Deletes feature with given probe_feature_ids
  Returntype : None
  Exceptions : Throws if not arg defines
  Caller     : Importer
  Status     : At Risk

=cut

sub delete_features{
	my ($self, $fids_ref) = @_;
	
	if(! @$fids_ref){
	  throw('Need to pass a ref to an array of feature ids');
	}
	
	my $cmd = 'DELETE from probe_feature WHERE probe_feature_id IN ('.join(',', @$fids_ref).')';
	$self->db->dbc->do($cmd);

	#This will fail anyway?
	#if($?){
	#  throw("SQL Command failed:\t$sql\n$@");
	#}

	return;
}

# All the results methods may be moved to a ResultAdaptor

=head2 fetch_results_by_channel_analysis

  Arg [1]    : int - Probe dbID
  Arg [2]    : int - Channel dbID
  Arg [1]    : string - Logic name of analysis
  Example    : my @results = @{$ofa->fetch_results_by_channel_analysis($op_id, $channel_id, 'RAW_VALUE')};
  Description: Gets all analysis results for probe on given channel
  Returntype : ARRAYREF
  Exceptions : warns if analysis is not valid in Channel context
  Caller     : ProbeFeature
  Status     : At Risk - rename fetch_results_by_probe_channel_analysis

=cut



sub fetch_results_by_channel_analysis{
	my ($self, $probe_id, $channel_id, $logic_name) = @_;

	throw("deprecated, use ResultSetAdaptor");

	
	#Will this always be RAW_VALUE?

	my %channel_metrics = (
						   RawValue => 1,
						  );


	if(! defined $probe_id || ! defined $channel_id) {
		throw("Need to define a valid probe and channel dbID");
	}
		

	my $analysis_clause = "";

	if($logic_name){
		if(exists $channel_metrics{$logic_name}){
			$analysis_clause = "AND a.logic_name = \"$logic_name\"";
		}else{
			warn("$logic_name is not a channel specific metric\nNo results returned\n");
			return;
		}
	}

	my $query = "SELECT r.score, a.logic_name from result r, analysis a where r.probe_id =\"$probe_id\" AND r.table_name=\"channel\" AND r.table_id=\"$channel_id\" AND r.analysis_id = a.analysis_id $analysis_clause";
	
	return $self->dbc->db_handle->selectall_arrayref($query);
}

=head2 fetch_results_by_Probe_Analysis_experimental_chip_ids

  Arg [1]    : Bio::EnsEMBL::Funcgen::Probe
  Arg [2]    : Bio:EnsEMBL::Analysis
  Arg [2]    : ARRAYREF - Bio::EnsEMBLExperimentalChip dbIDs
  Example    : my @results = @{$ofa->fetch_results_by_channel_analysis($probe, $analysis, \@ec_ids)};
  Description: Gets all analysis results for probe within a set of ExperimentalChips
  Returntype : ARRAYREF
  Exceptions : warns if analysis is not valid in ExperimentalChip context
  Caller     : ProbeFeature
  Status     : At Risk 

=cut

sub fetch_results_by_Probe_Analysis_experimental_chip_ids{
	my ($self, $probe, $anal, $chip_ids) = @_;

	throw("deprecated, use ResultSetAdaptor");


	my $logic_name = $anal->logic_name();
	my @cs = @$chip_ids;

	#warn "Fetching result for $probe_id, @cs, $logic_name";
	
	my $table_ids;
	my $table_name = "experimental_chip";

	my %chip_metrics = (
			    VSN_GLOG => 1,
			    SangerPCR =>1,
			   );

	#else no logic name or not a chip metric, then return channel and metric=?


	if(! defined $probe->dbID() || ! @$chip_ids) {
		throw("Need to define a valid probe and pass a listref of experimental chip dbIDs");
	}
		

	my $analysis_clause = ($logic_name) ? "AND a.logic_name = \"$logic_name\"" : "";

	if(! exists $chip_metrics{$logic_name}){
	  $table_name = "channel";
	  warn("Logic name($logic_name) is not a chip specific metric\nNo results returned\n");
	  
	  #build table ids from exp chip channel ids
	  #need to then sort out which channel is which in caller.

	  #need to enable raw data retrieval!!
	  return;
	}else{
	  $table_ids = join(", ", @$chip_ids);
	}


	#my $query = "SELECT r.score, r.table_id, a.logic_name from result r, analysis a where r.probe_id ='".$probe-dbID().
	#"' AND r.table_name=\"${table_name}\" AND r.table_id IN (${table_ids}) AND r.analysis_id = a.analysis_id $analysis_clause order by r.score";

	my $query = "SELECT r.score, r.table_id, a.logic_name from result r, analysis a where r.probe_id ='".$probe-dbID().
	  "' AND r.table_name=\"${table_name}\" AND r.table_id IN (${table_ids}) AND r.analysis_id = a.analysis_id $analysis_clause order by r.score";
	
	return $self->dbc->db_handle->selectall_arrayref($query);
}



sub _get_best_result{
  my ($self, $ofs, $analysis, $exp_chips) = @_;
  my ($median);

  throw('Deprecated, use EFGUtils');


  if(scalar(@$ofs) > 1){
    my @results = map $_->get_result_by_Analysis_ExperimentalChips($analysis, $exp_chips), @$ofs;
    @results = sort @results;

    my $count = scalar(@results);
    my $index = $count -1;
    #need to account for features/probes without results.  How would this happen?  Not all probes present in result file or score = NA?!
    #while(! $results[0]){
    #shift @results;
    #}
    if ($count == 1){
      $median =  $results[0];
    }
    elsif ($count % 2) { #odd number of scores
      $median = $results[($index+1)/2];
    }
    else { #even
      $median = ($results[($index)/2] + $results[($index/2)+1] ) / 2;
    }
  }else{
    $median = $ofs->[0]->get_result_by_Analysis_ExperimentalChips($analysis, $exp_chips);
  }

  return $median;
}




1;

