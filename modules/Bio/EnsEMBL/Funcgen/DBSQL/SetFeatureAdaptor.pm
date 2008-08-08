#
# Ensembl module for Bio::EnsEMBL::DBSQL::Funcgen::SetFeatureAdaptor
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::DBSQL::Funcgen::SetFeatureAdaptor - A base database adaptor for SetFeature adaptors.
storing SetFeature objects.

=head1 SYNOPSIS

my $afa = $db->get_SetFeatureAdaptor();

my $features = $afa->fetch_all_by_Slice($slice);


=head1 DESCRIPTION

The SetFeatureAdaptor is a base adaptor for all SetFeature adaptors.
e.g. SetFeature, AnnotatedFeature etc.  It provides common methods
across all feature types.

=head1 AUTHOR

This module was created by Nathan Johnson.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::DBSQL::SetFeatureAdaptor;

use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Funcgen::SetFeature;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor;

use vars qw(@ISA @EXPORT);
@ISA = qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor);
@EXPORT = (@{$DBI::EXPORT_TAGS{'sql_types'}});#required for child adaptor's store and _obj_from_sth methods





=head2 fetch_all_by_FeatureType_FeatureSets

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : Bio::EnsEMBL::FeatureType
  Arg [3]    : (optional) string - analysis logic name
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $features = $ofa->fetch_all_by_Slice_FeatureType($slice, $ft);
  Description: Retrieves a list of features on a given slice, specific for a given FeatureType.
  Returntype : Listref of Bio::EnsEMBL::SetFeature objects
  Exceptions : Throws if no FeatureType object provided
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_FeatureType_FeatureSets {
  my ($self, $ftype, @fsets) = @_;
	

  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType', $ftype);

  my @fs_ids;

  throw('Must provide a list of Bio::EnsEMBL::FeatureSet objects') if scalar(@fsets) == 0;

  foreach my $fset (@fsets) {
	$self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureSet', $fset);
	push (@fs_ids, $fset->dbID());
  }

  my $fs_ids = join(',', @fs_ids) if scalar(@fs_ids >1);

  my $constraint = $self->_main_table->[1].'.feature_set_id = fs.feature_set_id AND '.
	$self->_main_table->[1].'.feature_type_id='.$ftype->dbID.' AND '.$self->_main_table->[1].'.feature_set_id '; 
  $constraint .= (scalar(@fs_ids) >1) ? "IN ($fs_ids)" : '='.$fs_ids[0];

  return $self->generic_fetch($constraint);
}

=head2 fetch_all_by_Slice_FeatureType

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : Bio::EnsEMBL::FeatureType
  Arg [3]    : (optional) string - analysis logic name
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $features = $ofa->fetch_all_by_Slice_FeatureType($slice, $ft);
  Description: Retrieves a list of features on a given slice, specific for a given FeatureType.
  Returntype : Listref of Bio::EnsEMBL::SetFeature objects
  Exceptions : Throws if no FeatureType object provided
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_Slice_FeatureType {
  my ($self, $slice, $type, $logic_name) = @_;
	
  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType', $type);
  
  my $ft_id = $type->dbID();
  
  my $constraint = $self->_main_table->[1].".feature_set_id = fs.feature_set_id AND ".
	"fs.feature_type_id = '$ft_id'";

  $constraint = $self->_logic_name_to_constraint($constraint, $logic_name);

  return $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint);
}


=head2 fetch_all_by_FeatureSets

  Arg [2]    : Array of Bio::EnsEMBL::FeatureSet objects
  Example    : my $features = $set_feature_adaptor->fetch_all_by_FeatureSets(@fsets);
  Description: Retrieves a list of features specific for a given list of FeatureSets.
  Returntype : Listref of Bio::EnsEMBL::SetFeature objects
  Exceptions : Throws if list provided does not contain FeatureSets or if none provided
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_FeatureSets {
  my ($self, @fsets) = @_;
	
  my @fs_ids;

  throw('Must provide a list of Bio::EnsEMBL::FeatureSet objects') if scalar(@fsets) == 0;

  foreach my $fset (@fsets) {
	throw('Not a FeatureSet object') 
	  if ! (ref($fset) && $fset->isa("Bio::EnsEMBL::Funcgen::FeatureSet"));
	push (@fs_ids, $fset->dbID());
  }

  

  my $fs_ids = join(',', @fs_ids) if scalar(@fs_ids >1);
  my $constraint = $self->_main_table->[1].'.feature_set_id '; 
  $constraint .= (scalar(@fs_ids) >1) ? "IN ($fs_ids)" : '='.$fs_ids[0];

  #could have individual logic_names for each annotated feature here?
  #$constraint = $self->_logic_name_to_constraint($constraint, $logic_name);

  return $self->generic_fetch($constraint);
}




=head2 fetch_all_by_Slice_FeatureSets

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : Array of Bio::EnsEMBL::FeatureSet objects
  Example    : my $slice = $sa->fetch_by_region('chromosome', '1');
               my $features = $ofa->fetch_by_Slice_FeatureSets($slice, @fsets);
  Description: Retrieves a list of features on a given slice, specific for a given list of FeatureSets.
  Returntype : Listref of Bio::EnsEMBL::SetFeature objects
  Exceptions : Throws if list provided does not contain FeatureSets or if none provided
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_Slice_FeatureSets {
  my ($self, $slice, @fsets) = @_;
	
  my @fs_ids;


  throw('Must provide a list of Bio::EnsEMBL::FeatureSet objects') if scalar(@fsets) == 0;

  foreach my $fset (@fsets) {
	throw('Not a FeatureSet object') 
	  if ! ($fset && ref($fset) && $fset->isa("Bio::EnsEMBL::Funcgen::FeatureSet"));
	push (@fs_ids, $fset->dbID());
  }

  

  my $fs_ids = join(',', @fs_ids) if scalar(@fs_ids >1);
  my $constraint = $self->_main_table->[1].'.feature_set_id '; 
  $constraint .= (scalar(@fs_ids) >1) ? "IN ($fs_ids)" : '='.$fs_ids[0];

  #could have individual logic_names for each annotated feature here?
  #$constraint = $self->_logic_name_to_constraint($constraint, $logic_name);

  return $self->SUPER::fetch_all_by_Slice_constraint($slice, $constraint);
}



# Redefine BaseFeatureAdaptor method as analysis now abstracted to feature_set
# Given a logic name and an existing constraint this will
# add an analysis table constraint to the feature.  Note that if no
# analysis_id exists in the columns of the primary table then no
# constraint is added at all
# DO WE HAVE TO CALL THIS EXPLICITLY in this adaptor, or will BaseAdaptor use redefined method?


sub _logic_name_to_constraint {
  my $self = shift;
  my $constraint = shift;
  my $logic_name = shift;

  return $constraint if (!$logic_name);

  my $aa = $self->db->get_AnalysisAdaptor();
  my $an = $aa->fetch_by_logic_name($logic_name);

  if(! $an) {
	  #warn or throw?
	  warn("No analysis associated with logic_name $logic_name");
	  return undef;
  }

  my $an_id = $an->dbID();

  $constraint .= ' AND' if($constraint);
  $constraint .= " fs.analysis_id = $an_id";
  return $constraint;
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
	
  return $self->_main_table->[1].'.feature_set_id = fs.feature_set_id';
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
  Status     : At Risk

=cut


sub _final_clause {
  my $self = shift;
  return ' ORDER BY '.$self->_main_table->[1].'.seq_region_id, '.$self->_main_table->[1].'.seq_region_start';
}


#re-written for non-standard feature table

=head2 fetch_all_by_logic_name 

  Arg [1]     : string $logic_name the logic name of the analysis of features to obtain 
  Example     : $fs = $a->fetch_all_by_logic_name('foobar'); 
  Description : Returns an arrayref of features created from the database. only 
                features with an analysis of type $logic_name will be returned. 
  Returntype  : arrayref of Bio::EnsEMBL::SetFeatures
  Exceptions  : thrown if $logic_name not defined
  Caller      : General 
  Status      : At risk

=cut 

sub fetch_all_by_logic_name { 
  my $self = shift; 
  my $logic_name = shift || throw( "Need a logic_name" ); 

  my $constraint;
  my $an = $self->db->get_AnalysisAdaptor->fetch_by_logic_name($logic_name); 
  $constraint = ' '.$self->_main_table->[1].'.feature_set_id=fs.feature_set_id and fs.analysis_id = '.$an->dbID() if($an);

  return (defined $constraint) ? $self->generic_fetch($constraint) : undef;
} 


=head2 store_associated_feature_types

  Arg [1]     : string $logic_name the logic name of the analysis of features to obtain 
  Example     : $fs = $a->fetch_all_by_logic_name('foobar'); 
  Description : Returns an arrayref of features created from the database. only 
                features with an analysis of type $logic_name will be returned. 
  Returntype  : arrayref of Bio::EnsEMBL::SetFeatures
  Exceptions  : thrown if $logic_name not defined
  Caller      : General 
  Status      : At risk

=cut 

sub store_associated_feature_types { 
  my ($self, $set_feature) = @_;

  #Direct access to avoid lazy loading with an unstored SetFeature
  my $assoc_ftypes = $set_feature->{'associated_feature_types'};

  #Could be undef or empty
  return if ! defined $assoc_ftypes || scalar(@$assoc_ftypes) == 0;

  my $type = $set_feature->feature_set->type;
  my $feature_id = $set_feature->dbID;

  my $sql = 'INSERT into associated_feature_type(feature_id, feature_table, feature_type_id) values (?,?,?)';

  foreach my $ftype(@$assoc_ftypes){

	#We have already tested the class but not whether it is stored
	$self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType', $ftype);


	warn "$feature_id, $type, ".$ftype->dbID;

	my $sth = $self->prepare($sql);
	$sth->bind_param(1, $feature_id,  SQL_INTEGER);
	$sth->bind_param(2, $type,        SQL_VARCHAR);
	$sth->bind_param(3, $ftype->dbID, SQL_INTEGER);

	$sth->execute();
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
	
	return $self->_list_dbIDs($self->_main_table->[0]);
}

=head2 _main_table

  Example    : my $syn = $adaptor->_main_table->[1];
  Description: Convenience method to retrieve the main table or main table synonym for this adaptor
               Entirely dependent on ensembl convention of always having main table as first element
               of tables array.
  Returntype : Arrayref
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub _main_table{
  my $self = shift;

  #Need to do this to put it in list context to avoid just returning the last value
  my @tables = $self->_tables();
  return $tables[0];
}




1;
