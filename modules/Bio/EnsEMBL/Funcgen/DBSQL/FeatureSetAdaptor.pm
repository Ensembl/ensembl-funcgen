#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::FeatureSetAdaptor
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::Funcgen::DBSQL::FeatureSetAdaptor - A database adaptor for fetching and
storing Funcgen feature sets.

=head1 SYNOPSIS

my $fs_adaptor = $db->get_FeatureSetAdaptor();

my @fsets = $fs_adaptor->fetch_all_by_Experiment($exp);
my @displayable_fsets = @{$fs_adaptor->fetch_all_displayable()};

=head1 DESCRIPTION

The FeatureSetAdaptor is a database adaptor for storing and retrieving
Funcgen feature set.  

=head1 AUTHOR

This module was created by Nathan Johnson.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::DBSQL::FeatureSetAdaptor;

use Bio::EnsEMBL::Utils::Exception qw( warning );
use Bio::EnsEMBL::Funcgen::FeatureSet;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;

use vars qw(@ISA);


#May need to our this?
@ISA = qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor);



=head2 fetch_all_by_FeatureType

  Arg [1]    : Bio::EnsEMBL::Funcgen::FeatureType
  Arg [2]    : (optional) string - status e.g. 'DISPLAYABLE'
  Example    : my @fsets = $fs_adaptopr->fetch_all_by_FeatureType($type);
  Description: Retrieves FeatureSet objects from the database based on feature_type id.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::FeatureSet objects
  Exceptions : Throws if arg is not a valid FeatureType
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_FeatureType {
    my $self = shift;
    my $ftype = shift;
	my $status = shift;
    
	if(! ($ftype && $ftype->isa("Bio::EnsEMBL::Funcgen::FeatureType"))){
		throw("Must provide a valid Bio::EnsEMBL::Funcgen::FeatureType object");
	}
	
	my $sql = "fs.feature_type_id = '".$ftype->dbID()."'";
	$sql = $self->status_to_constraint($self, $constraint, $status) if $status;
	   

    return $self->generic_fetch("fs.feature_type_id = '".$ftype->dbID()."'");
	
}



=head2 fetch_all_by_CellType

  Arg [1]    : Bio::EnsEMBL::Funcgen::CelLType
  Arg [2]    : (optional) string - status e.g. 'DISPLAYABLE'
  Example    : my @displayable_fsets = @{$fset_adaptor->fetch_all_by_CellType($ctype, 'DISPLAYABLE')};
  Description: Fetch all FeatureSet with a given CellType.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::FeatureSet objects
  Exceptions : Throws if no valid CellType passed
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_CellType {
	my ($self, $ctype, $status) = @_;
	
	throw("Must provide a valid Bio::EnsEMBL::Funcgen::CellType") if (! ($ctype && $ctype->isa("Bio::EnsEMBL::Funcgen::CellType")));
	
	my $constraint = "fs.cell_type_id='".$ctype->dbID()."'";
	$constraint = $self->status_to_constraint($self, $constraint, $status) if $status;
	
	return $self->generic_fetch($constraint);
}

=head2 fetch_attributes

  Arg [1]    : Bio::EnsEMBL::Funcgen::OligoArray - array to fetch attributes for
  Example    : None
  Description: This function is solely intended to lazy load attributes into
               empty OligoArray objects. You should not need to call this.
  Returntype : None
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::OligoArray getters
  Status     : Medium Risk

=cut

sub fetch_attributes {
    my $self = shift;
    my $array = shift;

    my $tmp_array = $self->fetch_by_dbID( $array->dbID() );
    %$array = %$tmp_array;
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
	
	return ['feature_set', 'fs'];
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
	
	return qw( fs.feature_set_id fs.feature_type_id ft.analysis_id ft.feature_type_id);
}

=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates OligoArray objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::FeatureSet objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
	my ($self, $sth) = @_;
	
	my (@fsets, $fset, $analysis, %analysis_hash, $feature_type, $cell_type);
	my ($feature_set_id, $ftype_id, $analysis_id, $ctype_id, %ftype_hash);
	
	my $ft_adaptor = $self->db->get_FeatureTypeAdaptor();
	my $anal_adaptor = $self->db->get_AnalysisAdaptor();
	my $ct_adaptor = $self->db->get_CellTypeAdaptor();
	$ctype_hash{'NULL'} = undef;

	$sth->bind_columns();
	
	while ( $sth->fetch($feature_set_id, $ftype_id, $analysis_id, $ctype_id)) {

		$ctype_id ||= 'NULL';
		
		# Get the analysis object
		$analysis_hash{$analysis_id} = $anal_adaptor->fetch_by_dbID($analysis_id) if(! exists $analysis_hash{$analysis_id});

		# Get the feature type object
		$ftype_hash{$ftype_id} = $anal_adaptor->fetch_by_dbID($analysis_id) if(! exists $analysis_hash{$analysis_id});
		
		# Get the cell_type object
		$ctype_hash{$ctype_id} = $anal_adaptor->fetch_by_dbID($analysis_id) if(! exists $analysis_hash{$analysis_id});



		$fset = Bio::EnsEMBL::Funcgen::FeatureSet->new(
							       -dbID         => $feature_set_id,
							       -adaptor      => $self,
							       -feature_type => $ftype_hash{$ftype_id},
							       -analysis     => $analysis_hash{$analysis_id},
							       -cell_type    => $ctype_hash{$ctype_id},
							      );

		push @fsets, $fset;

	}

	return \@fsets;
}




=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::FeatureSet objects
  Example    : $oaa->store($fset1, $fset2, $fset3);
  Description: Stores FeatureSet objects in the database.
  Returntype : Listref of stored FeatureSet objects
  Exceptions : Throws if FeatureSet does not have a stored FeatureType
               Throws if invalid FeatureSet passed
  Caller     : General
  Status     : At Risk

=cut

sub store {
    my $self = shift;
    my @fsets = @_;

	#if(scalar(@fsets) == 0);

	my $sth = $self->prepare("INSERT INTO feature_set
                                 (feature_type_id, analysis_id, cell_type_id)
                                 VALUES (?, ?, ?)");


    foreach my $fset (@fsets) {
		throw('Can only store FeatureSet objects, skipping $fset')	if ( ! $fset->isa('Bio::EnsEMBL::Funcgen::FeatureSet'));

		if (!( $fset->dbID() && $fset->adaptor() == $self )){#use is_stored?

			
			#my $s_fset = $self->fetch_by_unique_and_experiment_id($ec->unique_id(), $ec->experiment_id());
			#throw("ExperimentalChip already exists in the database with dbID:".$s_ec->dbID().
			#	  "\nTo reuse/update this ExperimentalChip you must retrieve it using the ExperimentalChipAdaptor".
			#	  "\nMaybe you want to use the -recover option?") if $s_ec;
		

			 
			throw("FeatureSet must have a stored FeatureType") if (! $fset->feature_type->is_stored($self->db()));
			 
			my $ctype_id = (defined $fset->cell_type()) ? $fset->cell_type->dbID() : undef;

			$sth->bind_param(1, $fset->feature_type->dbID(), SQL_INTEGER);
			$sth->bind_param(2, $fset->analysis->dbID(),     SQL_INTEGER);
			$sth->bind_param(3, $ctype_id,                   SQL_INTEGER);
		
			$sth->execute();
			$fset->dbID($sth->{'mysql_insertid'});
			$fset->adaptor($self);
		}else{
			#assume we want to update the states
			warn('You may want to use $fset->adaptor->store_states($fset)');
			$self->store_states($fset);
		}
	}
	return \@fsets;
}

=head2 list_dbIDs

  Args       : None
  Example    : my @array_ids = @{$oaa->list_dbIDs()};
  Description: Gets an array of internal IDs for all OligoArray objects in the
               current database.
  Returntype : List of ints
  Exceptions : None
  Caller     : ?
  Status     : Medium Risk

=cut

sub list_dbIDs {
    my ($self) = @_;
	
    return $self->_list_dbIDs('feature_set');
}



1;

