#
# Ensembl module for Bio::EnsEMBL::DBSQL::Funcgen::DataSetAdaptor
#

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.


=head1 NAME

Bio::EnsEMBL::DBSQL::Funcgen::DataSetAdaptor - A database adaptor for fetching and
storing DataSet objects.

=head1 SYNOPSIS

my $dset_adaptor = $db->get_DataSetAdaptor();

my $dset = $dset_adaptor->fetch_by_dbID(1);
my @displayable_dsets = $dset_adaptor->fetch_all_displayable();

=head1 DESCRIPTION

The DataSetAdaptor is a database adaptor for storing and retrieving
DataSet objects.

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::DataSetAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Funcgen::DataSet;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;#sql_types bareword import

use base qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor);

=head2 fetch_by_name

  Arg [1]    : string - name of DataSet
  Arg [2]    : (optional) string - status e.g. 'DISPLAYABLE'
  Example    : my $dset = $dset_adaptor->fetch_by_name('RegulatoryFeatures:MultiCell');
  Description: Fetch DataSet with a given name
  Returntype : Bio::EnsEMBL::Funcgen::DataSet
  Exceptions : Throws if no name passed
  Caller     : General
  Status     : At Risk

=cut

sub fetch_by_name {
  my ($self, $name, $status) = @_;

  throw("Must provide a name argument") if (! defined $name);

  my $sql = "ds.name='${name}'";

  if ($status) {
    my $constraint = $self->status_to_constraint($status) if $status;
    $sql = (defined $constraint) ? $sql." ".$constraint : undef;
  }

  return (defined $sql) ? $self->generic_fetch($sql)->[0] : [];

}


=head2 fetch_all_by_supporting_set_type

  Arg [1]    : string - type of supporting_sets i.e. result or feature
  Arg [2]    : (optional) string - status e.g. 'DISPLAYABLE'
  Example    : my $dsets = $dset_adaptor->fetch_all_by_supporting_set('feature');
  Description: Fetch all DataSets whose pre-processed data consists of a particular set type
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::DataSet objects
  Exceptions : Throws if no supporting_set_type passed
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_supporting_set_type {
  my ($self, $type, $status) = @_;

  throw("Must provide a supporting_set type argument") if (! defined $type);

  my $sql = "ss.type='".$type."'";

  if ($status) {
    my $constraint = $self->status_to_constraint($status) if $status;
    $sql = (defined $constraint) ? $sql." ".$constraint : undef;
  }

  return (defined $sql) ? $self->generic_fetch($sql) : [];

}




=head2 fetch_by_product_FeatureSet

  Arg [1]    : Bio::EnsEMBL::Funcgen::FeatureSet
  Example    : my @dsets = $fs_adaptopr->fetch_by_product_FeatureSet($fset);
  Description: Retrieves DataSet objects from the database based on the FeatureSet.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::DataSet objects
  Exceptions : Throws if arg is not a valid FeatureSet
  Caller     : General
  Status     : At Risk

=cut

sub fetch_by_product_FeatureSet {
  my ($self, $fset) = @_;

  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureSet', $fset);
  $self->bind_param_generic_fetch($fset->dbID, SQL_INTEGER);

  return $self->generic_fetch('ds.feature_set_id = ?')->[0];
}


=head2 fetch_all_by_supporting_set

  Arg [1]    : Bio::EnsEMBL::Funcgen::Result|FeatureSet
  Example    : my @dsets = $fs_adaptopr->fetch_all_by_supporting_set($rset);
  Description: Retrieves DataSet objects from the database based on the
               given supporting Result or FeatureSet.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::DataSet objects
  Exceptions : Throws if arg is not a valid Result|FeatureSet
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_supporting_set {
  my ($self, $set) = @_;

  if (!
      (ref($set)
      &&
		  ($set->isa("Bio::EnsEMBL::Funcgen::ResultSet")  ||
		   $set->isa("Bio::EnsEMBL::Funcgen::FeatureSet") ||
		   $set->isa("Bio::EnsEMBL::Funcgen::InputSet")
      )
		  &&
     $set->dbID())) {
	throw("Must provide a valid stored Bio::EnsEMBL::Funcgen::ResultSet, FeatureSet or InputSet object");
  }

  #self join here to make sure we get all linked result_sets
  my $sql = ' ds.data_set_id IN (SELECT data_set_id from supporting_set where type="'.$set->set_type.'" and supporting_set_id='.$set->dbID().')';

#  warn "SQL $sql";



  return $self->generic_fetch($sql);
}


=head2 fetch_all_by_feature_type_class

  Arg [1]    : string - class of associated FeatureSet FeatureType
  Arg [2]    : optional: string - status e.g. 'DISPLAYABLE'
  Example    : my @dsets = @{$ds_adaptopr->fetch_all_by_feature_type_class('HISTONE')};
  Description: Retrieves DataSet objects from the database based on the FeatureSet FeatureType class.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::DataSet objects
  Exceptions : Throws if no class arg defined
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_feature_type_class {
  my ($self, $class, $status) = @_;

  throw ('Must provide a FeatureType class to retrieve DataSets') if ! defined $class;

  my ($constraint, @dsets);

  if ($status) {
	$constraint = $self->status_to_constraint($status) if $status;
	return [] if ! defined $constraint;
  }


  #This is fetching all feature sets!
  #we need to left join this?
  #we can't do it for class
  #but we can do it for product feature_set type

  foreach my $dset (@{$self->generic_fetch($constraint)}) {
	#uc both here to avoid case sensitivities
	push @dsets, $dset if uc($dset->product_FeatureSet->feature_type->class()) eq uc($class);
  }

  return \@dsets;
}


=head2 fetch_all_displayable_by_feature_type_class

  Arg [1]    : string - class of associated FeatureSet FeatureType
  Example    : my @dsets = @{$ds_adaptopr->fetch_all_displayable_by_feature_type_class('HISTONE')};
  Description: Wrapper method, retrieves all displayable DataSets with given FeatureSet FeatureType class
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::DataSet objects
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_displayable_by_feature_type_class {
  my ($self, $class) = @_;
  return $self->fetch_all_by_feature_type_class($class, 'DISPLAYABLE');
}


=head2 _true_tables

  Args       : None
  Example    : None
  Description: Returns the names and aliases of the tables to use for queries.
  Returntype : List of listrefs of strings
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _true_tables {
  return ([ 'data_set', 'ds' ], [ 'supporting_set', 'ss']);
}


=head2 _columns

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns a list of columns to use for queries.
  Returntype : List
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _columns {
  return qw(
			ds.data_set_id	    ds.feature_set_id
			ds.name             ss.type
			ss.supporting_set_id
		   );
}


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
  return (['supporting_set', 'ds.data_set_id = ss.data_set_id']);
}


=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates Array objects from an executed DBI statement
			   handle.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::DataSet objects
  Exceptions : Throws if supporting set cannot be fetched
  Caller     : Internal
  Status     : Stable

=cut

sub _objs_from_sth {
  my ($self, $sth) = @_;
  my (@data_sets, @supporting_sets, $data_set, $dbID, $set_id);
  my ($fset_id, $fset, $set, $name, $ss_type, $ss_id, %params);

  my %set_adaptors = (
	  feature      => $self->db->get_FeatureSetAdaptor,
	  result       => $self->db->get_ResultSetAdaptor,
    input        => $self->db->get_InputSetAdaptor,
	 );
  $sth->bind_columns(\$dbID, \$fset_id, \$name, \$ss_type, \$ss_id);

  while ( $sth->fetch() ) {

    if ((! %params) || ($params{-DBID} != $dbID)) {

	  if (%params) {
	   	push @data_sets, Bio::EnsEMBL::Funcgen::DataSet->new
		                    (%params, '-SUPPORTING_SETS', \@supporting_sets);
		  #do not set to empty array as this will cause failure of check in DataSet->new
		  undef @supporting_sets;
	  }

	  %params =(
	    -DBID                => $dbID,
	    -NAME                => $name,
	    -FEATURE_SET         => $set_adaptors{'feature'}->fetch_by_dbID($fset_id),
	    -ADAPTOR             => $self
     );
	}

	if ($ss_id) {
	  my $sset = $set_adaptors{$ss_type}->fetch_by_dbID($ss_id);

	  if (! $sset) {
		  throw('Could not find '.$name." supporting $ss_type set with dbID $ss_id");
	  }

	  push @supporting_sets, $sset;
	}
  }

  #handle last set
  if (%params) {
	push @data_sets, Bio::EnsEMBL::Funcgen::DataSet->new
		                   (%params, '-SUPPORTING_SETS', \@supporting_sets);
  }

  return \@data_sets;
}


=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::DataSet objects
  Example    : my @dsets = @{$dsa->store(@dsets)};
  Description: Stores given DataSet objects in the database. Sets dbID and adaptor
               on the objects that it stores.
  Returntype : ARRAYREF of stored DataSet objects
  Exceptions : Throws if no DataSet objects passed
               Throws if DataSet object has already been stored
               Throws if any supporting sets have not been stored
  Caller     : General
  Status     : At Risk

=cut

sub store{
  my ($self, @dsets) = @_;
  throw('Must pass a list of DataSet objects to store') if(! @dsets || $#dsets < 0);

  my $sth = $self->prepare("INSERT INTO data_set (feature_set_id, name) VALUES (?, ?)");
  my $sth2 = $self->prepare
    ("INSERT INTO supporting_set (data_set_id, supporting_set_id, type) VALUES (?, ?, ?)");
  my ($fset_id);
  my $db = $self->db();

  foreach my $dset (@dsets) {

    if( ! (defined $dset &&
          (ref($dset) eq 'Bio::EnsEMBL::Funcgen::DataSet'))){
      throw('Must pass a DataSet object to store')
    }

    if ( $dset->is_stored($db) ) {
      throw('DataSet [' . $dset->dbID() . '] is already stored in the database');
    }

    $sth->bind_param(1, $dset->product_FeatureSet->dbID, SQL_INTEGER);
    $sth->bind_param(2, $dset->name(),                   SQL_VARCHAR);
    $sth->execute();
    $dset->dbID( $self->last_insert_id );
    $dset->adaptor($self);

    foreach my $sset (@{$dset->get_supporting_sets()}) {

      if(! $sset->is_stored($db)){
        #Set check already done in _set_Sets_and_types
        throw('All supporting Feature and ResultSets must be stored previously');
      }

      $sth2->bind_param(1, $dset->dbID(),                SQL_INTEGER);
      $sth2->bind_param(2, $sset->dbID(),                SQL_INTEGER);
      $sth2->bind_param(3, $sset->set_type(),            SQL_VARCHAR); #enum feature/result/experimental
      $sth2->execute();
#       $self->store_states($dset);
    }
  }

  return \@dsets
}

1;

