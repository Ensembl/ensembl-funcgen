#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::EpigenomeAdaptor
#

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Funcgen::DBSQL::EpigenomeAdaptor - A database adaptor for fetching and
storing Funcgen Epigenome objects.

=head1 SYNOPSIS

my $eg_adaptor = $efgdba->get_EpigenomeAdaptor();

my $cell_type = $eg_adaptor->fetch_by_name("HeLa-S3");


=head1 DESCRIPTION

The EpigenomeAdaptor is a database adaptor for storing and retrieving
Funcgen Epigenome objects.

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::EpigenomeAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( warning throw );
use Bio::EnsEMBL::Funcgen::Epigenome;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;#sql_types barewords import

use base qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor);

#todo revert this to core BaseAdaptor and change _true_tables to _tables
# as we don't use query composition or other funcgen BaseAdaptor methods?

=head2 fetch_by_name

  Arg [1]    : string - name of Epigenome
  Arg [1]    : optional string - class of Epigenome
  Example    : my $eg = $eg_adaptor->fetch_by_name('HeLa');
  Description: Retrieves Epigenome objects by name.
  Returntype : Bio::EnsEMBL::Funcgen::Epigenome object
  Exceptions : Throws no name given
  Caller     : General
  Status     : At risk

=cut

sub fetch_by_name{
  my ($self, $name) = @_;
  throw("Must specify an Epigenome name") if ! defined $name;
  $self->bind_param_generic_fetch($name, SQL_VARCHAR);

  #name is unique so we should only have one
  return $self->generic_fetch('eg.name = ?')->[0];
}

sub fetch_by_display_label {
  my ($self, $name) = @_;
  throw("Must specify the display label of an Epigenome") if ! defined $name;
  $self->bind_param_generic_fetch($name, SQL_VARCHAR);

  #name is unique so we should only have one
  return $self->generic_fetch('eg.display_label = ?')->[0];
}

sub fetch_by_dbID {
  my ($self, $dbID) = @_;
  
  if (defined $self->{_cache}->{$dbID}) {
    return $self->{_cache}->{$dbID};
  }
  
  my $epigenome = $self->SUPER::fetch_by_dbID($dbID);
  $self->{_cache}->{$epigenome->dbID} = $epigenome;
  
  return $self->{_cache}->{$dbID};
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
  return (['epigenome', 'eg']);
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
  return qw( 
     eg.epigenome_id 
     eg.name 
     eg.display_label
     eg.description 
     eg.gender 
     eg.ontology_accession 
     eg.tissue 
     eg.production_name
  );
  #type/class = enum cell, cell line, tissue
}

=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates Channel objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::Epigenome objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
    my ( $self, $sth ) = @_;

    my (@result,             $eg_id,  $name,
        $dlabel,             $desc,   $gender,
        $ontology_accession, $tissue, $production_name
    );

    $sth->bind_columns( \$eg_id, \$name, \$dlabel, \$desc, \$gender,
        \$ontology_accession, \$tissue, \$production_name );

    while ( $sth->fetch() ) {
        my $epigenome = Bio::EnsEMBL::Funcgen::Epigenome->new(
            -dbID               => $eg_id,
            -NAME               => $name,
            -DISPLAY_LABEL      => $dlabel,
            -DESCRIPTION        => $desc,
            -GENDER             => $gender,
            -ONTOLOGY_ACCESSION => $ontology_accession,
            -TISSUE             => $tissue,
            -ADAPTOR            => $self,
            -PRODUCTION_NAME    => $production_name,
        );

        push @result, $epigenome;

    }
    return \@result;
}




=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::Epigenome objects
  Example    : $chan_a->store($c1, $c2, $c3);
  Description: Stores Epigenome objects in the database.
  Returntype : None
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub store {
  my $self = shift;
  my @args = @_;


  my $sth = $self->prepare("
			INSERT INTO epigenome
			(name, display_label, description, gender, ontology_accession, tissue, production_name)
			VALUES (?, ?, ?, ?, ?, ?, ?)");



  foreach my $eg (@args) {
	  if ( ! $eg->isa('Bio::EnsEMBL::Funcgen::Epigenome') ) {
		  warning('Can only store Epigenome objects, skipping $eg');
		  next;
	  }

	  if ( $eg->dbID() && $eg->adaptor() == $self ){
		  warn("Skipping previously stored Epigenome dbID:".$eg->dbID().")");
		  next;
	  }


    $sth->bind_param( 1, $eg->name,               SQL_VARCHAR );
    $sth->bind_param( 2, $eg->display_label,      SQL_VARCHAR );
    $sth->bind_param( 3, $eg->description,        SQL_VARCHAR );
    $sth->bind_param( 4, $eg->gender,             SQL_VARCHAR );
    $sth->bind_param( 5, $eg->ontology_accession, SQL_VARCHAR );
    $sth->bind_param( 6, $eg->tissue,             SQL_VARCHAR );
    $sth->bind_param( 7, $eg->production_name,    SQL_VARCHAR );

	  
	  $sth->execute();
	  $eg->dbID($self->last_insert_id);
	  $eg->adaptor($self);
	}

  return \@args;
}




1;

