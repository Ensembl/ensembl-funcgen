#
# Ensembl module for Bio::EnsEMBL::DBSQL::Funcgen::PredictedFeatureAdaptor
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::DBSQL::Funcgen::PredictedFeatureAdaptor - A database adaptor for fetching and
storing PredictedFeature objects.

=head1 SYNOPSIS

my $pfa = $db->get_PredictedFeatureAdaptor();

my $features = $pfa->fetch_all_by_Slice($slice);
$features = $pfa->fetch_all_by_Slice_Target($slice, $target);

=head1 DESCRIPTION

The PredictedFeatureAdaptor is a database adaptor for storing and retrieving
PredictedFeature objects.

=head1 AUTHOR

This module was created by Nathan Johnson.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::DBSQL::PredictedFeatureAdaptor;

use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Funcgen::PredictedFeature;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor);

#fetch_all_by_Slice_Experiment
#fetch_all_by_Slice_Analysis handled by BaseFeatureAdaptor fetch_all_by_Slice_score, which also takes logic name
#make wrapper method for useability?

#re-write the methods to take names aswell as objects? 


=head2 fetch_all_by_Slice_FeatureType

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : Bio::EnsEMBL::FeatureType
  Arg [3]    : (optional) string - logic name
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $features = $ofa->fetch_by_Slice_Target($slice, $target);
  Description: Retrieves a list of features on a given slice, specific for a given target.
  Returntype : Listref of Bio::EnsEMBL::PredictedFeature objects
  Exceptions : Throws if no Target object provided
  Caller     : General
  Status     : Medium Risk

=cut

sub fetch_all_by_Slice_FeatureType {
  my ($self, $slice, $type, $logic_name) = @_;
	
  throw('Need type as parameter') if ! $type->isa("Bio::EnsEMBL::FeatureType");
  my $ft_id = $type->dbID();
  my $constraint = qq( pf.feature_type_id = '$ft_id' );
  
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
			[ 'predicted_feature', 'pf' ],
	       );
			#[ 'experiment_prediction', 'ep'],
			
			#target?
		#   );
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
	    pf.predicted_feature_id  pf.seq_region_id
	    pf.seq_region_start      pf.seq_region_end
	    pf.seq_region_strand     pf.coord_system_id
	    pf.feature_type_id       pf.display_label
	    pf.analysis_id	     pf.score
	   );
}


#change this to experiment_prediction query?
#do we need this?

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
	
	#return 'pf.predicted_feature_id = ep.predicted_feature_id';

	return;
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
	return ' ORDER BY pf.seq_region_id, pf.seq_region_start, pf.predicted_feature_id';
}


=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates PredictedFeature objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::PredictedFeature objects
  Exceptions : None
  Caller     : Internal
  Status     : Medium Risk

=cut

sub _objs_from_sth {
	my ($self, $sth, $mapper, $dest_slice) = @_;

	#For EFG this has to use a dest_slice from core/dnaDB whether specified or not.
	#So if it not defined then we need to generate one derived from the species_name and schema_build of the feature we're retrieving.

	# This code is ugly because caching is used to improve speed

	#my $sa = $self->db->get_SliceAdaptor();
	
	my ($sa, $old_cs_id);
	$sa = $dest_slice->adaptor->db->get_SliceAdaptor() if($dest_slice);#don't really need this if we're using DNADBSliceAdaptor?

	#Some of this in now probably overkill as we'll always be using the DNADB as the slice DB
	#Hence it should always be on the same coord system

	my $aa = $self->db->get_AnalysisAdaptor();
	my @features;
	my (%analysis_hash, %slice_hash, %sr_name_hash, %sr_cs_hash);

	my (
	    $predicted_feature_id,  $seq_region_id,
	    $seq_region_start,  $seq_region_end,
	    $seq_region_strand, $cs_id,
	    $ft_id,        $display_label,    
	    $analysis_id,       $score,
	);
	$sth->bind_columns(
			   \$predicted_feature_id,  \$seq_region_id,
			   \$seq_region_start,  \$seq_region_end,
			   \$seq_region_strand, \$cs_id,
			   \$ft_id,        \$display_label,
			   \$analysis_id,		\$score,
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



	    #Need to build a slice adaptor cache here?
	    #Would only ever want to do this if we enable mapping between assemblies??
	    #Or if we supported the mapping between cs systems for a given schema_build, which would have to be handled by the core api
	    
	    
	    if($old_cs_id && ($old_cs_id != $cs_id)){
	      throw("More than one coord_system for feature query, need to implement SliceAdaptor hash?");
	    }
	    
	    $old_cs_id = $cs_id;
	    
	    
	    #Need to make sure we are restricting calls to Experiment and channel(i.e. the same coord_system_id)
	    
	    $sa ||= $self->db->get_SliceAdaptor($cs_id);
	    
	    
	    
	    # This assumes that features come out sorted by ID
	    next if ($last_feature_id == $predicted_feature_id);
	    $last_feature_id = $predicted_feature_id;
	    
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
	    

	
	    push @features, $self->_new_fast( {
					       'start'         => $seq_region_start,
					       'end'           => $seq_region_end,
					       'strand'        => $seq_region_strand,
					       'slice'         => $slice,
					       'analysis'      => $analysis,
					       'adaptor'       => $self,
					       'dbID'          => $predicted_feature_id,
					       'score'         => $score,
					       'display_label' => $display_label,
					       'feature_type_id'   => $ft_id
					      } );
	  }
	
	return \@features;
      }

=head2 _new_fast

  Args       : Hashref to be passed to PredictedFeature->new_fast()
  Example    : None
  Description: Construct a PredictedFeature object using quick and dirty new_fast.
  Returntype : Bio::EnsEMBL::Funcgen::PredictedFeature
  Exceptions : None
  Caller     : _objs_from_sth
  Status     : Medium Risk

=cut

sub _new_fast {
	my $self = shift;
	
	my $hash_ref = shift;
	return Bio::EnsEMBL::Funcgen::PredictedFeature->new_fast($hash_ref);
}

=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::PredictedFeature objects
  Example    : $ofa->store(@features);
  Description: Stores given PredictedFeature objects in the database. Should only
               be called once per feature because no checks are made for
			   duplicates. Sets dbID and adaptor on the objects that it stores.
  Returntype : None
  Exceptions : Throws if a list of PredictedFeature objects is not provided or if
               an analysis is not attached to any of the objects
  Caller     : General
  Status     : Medium Risk

=cut

sub store{
  my ($self, @pfs) = @_;
  
  if (scalar(@pfs) == 0) {
    throw('Must call store with a list of PredictedFeature objects');
  }

  my $sth = $self->prepare("
		INSERT INTO predicted_feature (
			seq_region_id,  seq_region_start,
			seq_region_end, seq_region_strand,
                        coord_system_id, feature_type_id,
			display_label,  analysis_id,
			score
		) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
	");

  my $db = $self->db();
  my $analysis_adaptor = $db->get_AnalysisAdaptor();
  
 FEATURE: foreach my $pf (@pfs) {
    
    if( !ref $pf || !$pf->isa('Bio::EnsEMBL::Funcgen::PredictedFeature') ) {
      throw('Feature must be an PredictedFeature object');
    }
    
    if ( $pf->is_stored($db) ) {
      warning('PredictedFeature [' . $pf->dbID() . '] is already stored in the database');
      next FEATURE;
    }

    if ( !defined $pf->analysis() ) {
      throw('An analysis must be attached to the PredictedFeature objects to be stored.');
    }

    # Store the analysis if it has not been stored yet
    $analysis_adaptor->store( $pf->analysis()) if ( !$pf->analysis->is_stored($db) );
    #add ft here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    my $original = $pf;
    my $seq_region_id;
    ($pf, $seq_region_id) = $self->_pre_store($pf);
    
    $sth->bind_param(1, $seq_region_id,         SQL_INTEGER);
    $sth->bind_param(2, $pf->start(),           SQL_INTEGER);
    $sth->bind_param(3, $pf->end(),             SQL_INTEGER);
    $sth->bind_param(4, $pf->strand(),          SQL_TINYINT);
    $sth->bind_param(5, $pf->coord_system_id(), SQL_INTEGER);
    $sth->bind_param(6, $pf->feature_type_id(), SQL_INTEGER);
    $sth->bind_param(7, $pf->display_label(),   SQL_VARCHAR);
    $sth->bind_param(8, $pf->analysis->dbID(),  SQL_INTEGER);
    $sth->bind_param(9, $pf->score(),            SQL_DOUBLE);
    
    $sth->execute();

    $original->dbID( $sth->{'mysql_insertid'} );
    $original->adaptor($self);
  }

  return;
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
	
	return $self->_list_dbIDs('predicted_feature');
}







1;

