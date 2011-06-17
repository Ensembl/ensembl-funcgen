#
# Ensembl module for Bio::EnsEMBL::DBSQL::Funcgen::ExperimentalSetAdaptor
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::DBSQL::Funcgen::ExperimentalSetAdaptor - A database adaptor for fetching and
storing ExperimentalSet objects.  

=head1 SYNOPSIS

my $rset_adaptor = $db->get_ExperimentalSetAdaptor();

my @rsets = @{$rset_adaptor->fetch_all_ExperimentalSets_by_Experiment()};
my @displayable_rsets = @{$rset_adaptor->fetch_all_displayable_ExperimentalSets()};

#Other methods?
#by FeatureType, CellType all with displayable flag?


=head1 DESCRIPTION

The ExperimentalSetAdaptor is a database adaptor for storing and retrieving
ExperimentalSet objects.

=head1 AUTHOR

This module was created by Nathan Johnson.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::DBSQL::ExperimentalSetAdaptor;

use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Funcgen::ExperimentalSet;
use Bio::EnsEMBL::Funcgen::ResultFeature;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(mean median);

use vars qw(@ISA);


@ISA = qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor);



=head2 fetch_all_by_FeatureType

  Arg [1]    : Bio::EnsEMBL::Funcgen::FeatureType
  Example    : 
  Description: Retrieves a list of features on a given slice that are created
               by probes from the specified type of array.
  Returntype : Listref of Bio::EnsEMBL::OligoFeature objects
  Exceptions : Throws if no array type is provided
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_FeatureType {
  my ($self, $ftype) = @_;

  if( !(ref($ftype) && $ftype->isa("Bio::EnsEMBL::Funcgen::FeatureType") && $ftype->dbID())){
    throw("Need to pass a valid stored Bio::EnsEMBL::Funcgen::FeatureType");
  }
  
  my $constraint = "es.feature_type_id =".$ftype->dbID();
	
  return $self->generic_fetch($constraint);
}


=head2 fetch_all_by_CellType

  Arg [1]    : Bio::EnsEMBL::Funcgen::CellType
  Example    : 
  Description: 
  Returntype : Arrayref of Bio::EnsEMBL::Funcgen::ExperimentalSet objects
  Exceptions : Throws if no CellType is provided
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_CellType {
  my ($self, $ctype) = @_;

  if( !(ref($ctype) && $ctype->isa("Bio::EnsEMBL::Funcgen::CellType") && $ctype->dbID())){
    throw("Need to pass a valid stored Bio::EnsEMBL::Funcgen::CellType");
  }
	
  my $constraint = "es.cell_type_id =".$ctype->dbID();
	
  return $self->generic_fetch($constraint);
}
 

=head2 fetch_all_by_Experiment

  Arg [1]    : Bio::EnsEMBL::Funcgen::Experiment
  Example    : $exp_set = $eseta->fetch_by_Experiment($exp);
  Description: Retrieves a ExperimentalSet based on the given Experiment
  Returntype : Bio::EnsEMBL::Funcgen::ExperimentalSet
  Exceptions : Throws if no valid stored Experiment provided
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_Experiment {
  my ($self, $exp) = @_;

  if( ! ( ref($exp) &&
		  $exp->isa('Bio::EnsEMBL::Funcgen::Experiment') &&
		  $exp->dbID())){
	throw('Need to pass a valid stored Bio::EnsEMBL::Funcgen::Experiment');
  }
		
  return $self->generic_fetch('es.experiment_id = '.$exp->dbID());
}

=head2 fetch_by_name

  Arg [1]    : string - ExperimentalSet name
  Example    : $exp_set = $eseta->fetch_by_Experiment('exp_set_1');
  Description: Retrieves a ExperimentalSet based on the ExperimetnalSet name
  Returntype : Bio::EnsEMBL::Funcgen::ExperimentalSet
  Exceptions : Throws if no name provided
  Caller     : General
  Status     : At Risk

=cut

sub fetch_by_name {
  my ($self, $name) = @_;

  throw('Need to pass a name argument') if( ! defined $name);
		
  return $self->generic_fetch("es.name ='${name}'")->[0];
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
		  [ 'experimental_set',    'es' ],
		  [ 'experimental_subset', 'ess' ],
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

	return qw(
			  es.experimental_set_id  es.experiment_id
			  es.feature_type_id      es.cell_type_id
			  es.format               es.vendor
			  es.name     			  ess.name
			  ess.experimental_subset_id
		 );

	
}

#=head2 _default_where_clause
#
#  Args       : None
#  Example    : None
#  Description: PROTECTED implementation of superclass abstract method.
#               Returns an additional table joining constraint to use for
#			   queries.
#  Returntype : List of strings
#  Exceptions : None
#  Caller     : Internal
#  Status     : At Risk
#
#=cut

#sub _default_where_clause {
#  my $self = shift;

#  return 'es.experimental_set_id = ess.experimental_set_id';

#}

=head2 _left_join

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns an additional table joining constraint to use for
			   queries.
  Returntype : List
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _left_join {
  my $self = shift;
	
  return (['experimental_subset', 'es.experimental_set_id = ess.experimental_set_id']);
}





=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates Array objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::Experiment objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
  my ($self, $sth) = @_;
  
  my ($dbid, $exp_id, $ftype_id, $ctype_id, $format, $vendor, $name, $ess_name, $ess_id);
  my ($eset, @esets, $ftype, $ctype);
  my $ft_adaptor = $self->db->get_FeatureTypeAdaptor();
  my $ct_adaptor = $self->db->get_CellTypeAdaptor();
  my $exp_adaptor = $self->db->get_ExperimentAdaptor();
  $sth->bind_columns(\$dbid, \$exp_id, \$ftype_id, \$ctype_id, \$format, \$vendor, \$name, \$ess_name, \$ess_id);
  
  #this fails if we delete entries from the joined tables
  #causes problems if we then try and store an rs which is already stored

  while ( $sth->fetch() ) {

    if(! $eset || ($eset->dbID() != $dbid)){
      
      push @esets, $eset if $eset;
      $ftype = (defined $ftype_id) ? $ft_adaptor->fetch_by_dbID($ftype_id) : undef;
      $ctype = (defined $ctype_id) ? $ct_adaptor->fetch_by_dbID($ctype_id) : undef;

      $eset = Bio::EnsEMBL::Funcgen::ExperimentalSet->new(
														  -DBID         => $dbid,
														  -EXPERIMENT   => $exp_adaptor->fetch_by_dbID($exp_id),
														  -FORMAT       => $format,
														  -VENDOR       => $vendor,
														  -FEATURE_TYPE => $ftype,
														  -CELL_TYPE    => $ctype,
														  -ADAPTOR      => $self,
														  -NAME         => $name,
														 );
    }
    
    #This assumes logical association between chip from the same exp, confer in store method?????????????????
	
	
	#we're not controlling ctype and ftype during creating new ExperimentalSets to store.
	#we should change add_table_id to add_ExperimentalChip and check in that method
    if(defined $ess_name){
	  
	  $eset->add_new_subset($ess_name, Bio::EnsEMBL::Funcgen::ExperimentalSubset->new( -name    => $ess_name,
																				   -dbID    => $ess_id,
																				   -adaptor => $self,
																				   -experimental_set => $eset,
																				 ));
	  
	}
  }
  
  push @esets, $eset if $eset;
  
  return \@esets;
}



=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::ExperimentalSet objects
  Example    : $rsa->store(@esets);
  Description: Stores or updates previously stored ExperimentalSet objects in the database. 
  Returntype : None
  Exceptions : Throws if a List of ExperimentalSet objects is not provided or if
               an analysis is not attached to any of the objects
  Caller     : General
  Status     : At Risk

=cut

sub store{
  my ($self, @exp_sets) = @_;

  throw("Must provide a list of ExperimentalSet objects") if(scalar(@exp_sets == 0));
  
  
  
  my $sth = $self->prepare('INSERT INTO experimental_set (experiment_id, feature_type_id, 
                                                       cell_type_id, format, vendor, name) 
                                                       VALUES (?, ?, ?, ?, ?, ?)');
  
  my $db = $self->db();
  
  foreach my $set (@exp_sets) {
    
    if( ! ref $set || ! $set->isa('Bio::EnsEMBL::Funcgen::ExperimentalSet') ) {
      throw('Must be an ExperimentalSet object to store');
    }
    
        
    if ( $set->is_stored($db) ) {
      throw('ExperimentalSet [' . $set->dbID() . '] is already stored in the database\nExperimentalSetAdaptor does not yet accomodate updating ExperimentalSets');
      #would need to retrive stored result set and update table_ids
    }
   

	my $ct_id = (defined $set->cell_type()) ? $set->cell_type->dbID() : undef;
	my $ft_id = (defined $set->feature_type()) ? $set->feature_type->dbID() : undef;

    $sth->bind_param(1, $set->get_Experiment->dbID(),   SQL_INTEGER);
	$sth->bind_param(2, $ft_id,                         SQL_INTEGER);
	$sth->bind_param(3, $ct_id,                         SQL_INTEGER);
  	$sth->bind_param(4, $set->format,                   SQL_VARCHAR);
  	$sth->bind_param(5, $set->vendor,                   SQL_VARCHAR);
	$sth->bind_param(6, $set->name,                     SQL_VARCHAR);
	
    
    $sth->execute();
    
    $set->dbID( $sth->{'mysql_insertid'} );
    $set->adaptor($self);

    
    $self->store_ExperimentalSubsets($set->get_subsets()) if @{$set->get_subsets()};
  }
  
  return \@exp_sets;
}


=head2 store_ExperimentalSubsets

  Args       : Bio::EnsEMBL::Funcgen::ExperimentalSet 
  Example    : $esa->store_ExperimentalSubsets(\@e_subsets);
  Description: Convenience methods extracted from store to allow updating of ExperimentalSubset entries 
               during inline result processing which would otherwise be troublesome due to the need
               for an ExperimentalSet 
  Returntype : Bio::EnsEMBL::Funcgen::ExperimentalSet
  Exceptions : Throws if a stored ExperimentalSet object is not provided
               Throws if no ExperimentalSubsets present
  Caller     : General
  Status     : At Risk

=cut


sub store_ExperimentalSubsets{
  my ($self, $ssets) = @_;
  
  my $sth = $self->prepare("
		INSERT INTO experimental_subset (
			experimental_set_id, name
		) VALUES (?, ?)
	");

  throw('Must provide at least one ExperimentalSubset') if(! @$ssets);

  #Store and set all previously unstored table_ids
  foreach my $sset(@$ssets){
	
	#use is_stored here?
	if($sset->dbID()){
	  warn "Skipping ExperimentalSubset ".$sset->name()." - already stored in the DB";
	  next;
	}
	

	$sth->bind_param(1, $sset->experimental_set->dbID(), SQL_INTEGER);
	$sth->bind_param(2, $sset->name(),                   SQL_VARCHAR);
	$sth->execute();

	$sset->dbID($sth->{'mysql_insertid'});
	$sset->adaptor($self);


	#No need to set it as we're working on the hasref here, so should be updated in the class.
	#add directly to avoid name clash warnings
	#$exp_set->{'subsets'}{$sub_set_name} = Bio::EnsEMBL::Funcgen::ExperimentalSubset->new
	#  (
	#   -dbID    => $sth->{'mysql_insertid'},
	#   -name    => $sub_set_name,
	#   -adaptor => $self,
	#   #-experimental_set_id?
	#  );

  }
  
  #don't really need to return as we're passing the ref
  return $ssets;
}

=head2 list_dbIDs

  Args       : None
  Example    : my @sets_ids = @{$esa->list_dbIDs()};
  Description: Gets an array of internal IDs for all ExperimentalSet objects in
               the current database.
  Returntype : List of ints
  Exceptions : None
  Caller     : general
  Status     : stable

=cut

sub list_dbIDs {
	my $self = shift;
	
	return $self->_list_dbIDs('result_set');
}

1;

