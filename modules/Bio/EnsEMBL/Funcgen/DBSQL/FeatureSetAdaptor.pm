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

Bio::EnsEMBL::Funcgen::DBSQL::FeatureSetAdaptor - A database adaptor for fetching and
storing Funcgen feature sets.

=head1 SYNOPSIS

my $fs_adaptor = $db->get_FeatureSetAdaptor();

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::FeatureSetAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw deprecate );
use Bio::EnsEMBL::Funcgen::FeatureSet;
use Bio::EnsEMBL::Funcgen::DBSQL::SetAdaptor; #DBI sql_types import

use base qw(Bio::EnsEMBL::Funcgen::DBSQL::SetAdaptor);

=head2 fetch_by_name

  Arg [1]    : string - name of FeatureSet
  Arg [2]    : (optional) string - status e.g. 'DISPLAYABLE'
  Example    : my @fsets = @{$fset_adaptor->fetch_by_name('feature_set-1')};
  Description: Fetch all FeatureSets wit a given name
  Returntype : Bio::EnsEMBL::Funcgen::FeatureSet objects
  Exceptions : Throws if no name passed
  Caller     : General
  Status     : At Risk

=cut

sub fetch_by_name {
  my ($self, $name, $status) = @_;

  throw("Must provide a name argument") if (! defined $name);

  my $sql = "fs.name='".$name."'";

  if($status){
    my $constraint = $self->status_to_constraint($status) if $status;
    $sql = (defined $constraint) ? $sql." ".$constraint : undef;
  }

  return (defined $sql) ? $self->generic_fetch($sql)->[0] : [];

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
  return ([ 'feature_set', 'fs' ]);
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

	return qw( fs.feature_set_id fs.feature_type_id
			   fs.analysis_id fs.name fs.type
			   fs.description fs.display_label);
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

	my (@fsets, $fset, $analysis, %analysis_hash, $feature_type, $epigenome, $name, $type, $display_label, $desc);
	my ($feature_set_id, $ftype_id, $analysis_id, %ftype_hash, %epigenome_hash);

	my $ft_adaptor = $self->db->get_FeatureTypeAdaptor();
	my $anal_adaptor = $self->db->get_AnalysisAdaptor();
	$epigenome_hash{'NULL'} = undef;

	$sth->bind_columns(\$feature_set_id, \$ftype_id, \$analysis_id,
                       \$name, \$type, \$desc, \$display_label,);

	while ( $sth->fetch()) {

		# Get the analysis object
		$analysis_hash{$analysis_id} = $anal_adaptor->fetch_by_dbID($analysis_id) if(! exists $analysis_hash{$analysis_id});

		# Get the feature type object
		$ftype_hash{$ftype_id} = $ft_adaptor->fetch_by_dbID($ftype_id) if(! exists $ftype_hash{$ftype_id});

		#Use new_fast here and strip the prefixed -'s
		$fset = Bio::EnsEMBL::Funcgen::FeatureSet->new
		  (
		   -dbID          => $feature_set_id,
		   -adaptor       => $self,
		   -feature_type  => $ftype_hash{$ftype_id},
		   -analysis      => $analysis_hash{$analysis_id},
		   -name          => $name,
		   -feature_class => $type,
		   -display_label => $display_label,
		   -description   => $desc,
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
               Throws if not FeatureSets passed
               Warns if external_db_name not defined is type is external
               Throws if external_db is not present in the db
  Caller     : General
  Status     : At Risk

=cut

sub store {
    my $self = shift;
    my @fsets = @_;

	throw('Must supply a list of FeatureSets to store') if(scalar(@fsets) == 0);

	my $sth = $self->prepare
	  (
	   "INSERT INTO feature_set
        (feature_type_id, analysis_id, name, type, description, display_label)
        VALUES (?, ?, ?, ?, ?, ?)"
	  );

	my $db = $self->db;
	my ($sql, $edb_id, %edb_hash);

    foreach my $fset (@fsets) {
		throw('Can only store FeatureSet objects, skipping $fset')	if ( ! $fset->isa('Bio::EnsEMBL::Funcgen::FeatureSet'));


		if (! $fset->is_stored($db) ) {

		  # Check FeatureType and Analysis
		  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType', $fset->feature_type);
		  $self->db->is_stored_and_valid('Bio::EnsEMBL::Analysis', $fset->analysis);

		  $sth->bind_param(1, $fset->feature_type->dbID, SQL_INTEGER);
		  $sth->bind_param(2, $fset->analysis->dbID,     SQL_INTEGER);
		  $sth->bind_param(3, $fset->name,               SQL_VARCHAR);
		  $sth->bind_param(4, $fset->feature_class,      SQL_VARCHAR);
		  $sth->bind_param(5, $fset->description,        SQL_VARCHAR);
		  $sth->bind_param(6, $fset->display_label,      SQL_VARCHAR);

		  $sth->execute();
		  $fset->dbID($self->last_insert_id);
		  $fset->adaptor($self);
		}
		else{
			warn('FeatureSet '.$fset->name.'is already stored, updating status entries');
			$self->store_states($fset);
		}
	}
	return \@fsets;
}

1;



