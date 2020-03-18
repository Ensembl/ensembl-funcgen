#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::EpigenomeAdaptor
#

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(stringify destringify);

use base qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor);
use Bio::EnsEMBL::Utils::Exception qw( throw warning deprecate );

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

sub fetch_by_production_name {
  my ($self, $name) = @_;
  throw("Must specify an Epigenome name") if ! defined $name;
  $self->bind_param_generic_fetch($name, SQL_VARCHAR);

  # production_name is unique so we should only have one
  return $self->generic_fetch('eg.production_name = ?')->[0];
}

sub fetch_by_short_name {
  my ($self, $name) = @_;
  throw("Must specify the short name of an Epigenome") if ! defined $name;
  $self->bind_param_generic_fetch($name, SQL_VARCHAR);

  #short_name is unique so we should only have one
  return $self->generic_fetch('eg.short_name = ?')->[0];
}

sub fetch_by_dbID {
  my ($self, $dbID) = @_;
  
  if (defined $self->{_cache}->{$dbID}) {
    return $self->{_cache}->{$dbID};
  }
  
  my $epigenome = $self->SUPER::fetch_by_dbID($dbID);
  
  if (! defined $epigenome) {
    use Carp;
    confess("Can't fetch epigenome with epigenome_id $dbID!");
  }
  
  $self->{_cache}->{$epigenome->dbID} = $epigenome;
  
  return $epigenome;
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
     eg.short_name
     eg.description 
     eg.gender 
     eg.production_name
     eg.search_terms
     eg.full_name
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
        $short_name,             $desc,   $gender,
        $production_name, $search_terms, $full_name,
    );

    $sth->bind_columns( \$eg_id, \$name, \$short_name, \$desc, \$gender,
        \$production_name, \$search_terms, \$full_name);

    while ( $sth->fetch() ) {

        my $search_terms_ref;
        if (defined $search_terms){
          $search_terms_ref = destringify($search_terms);
        }

        my $epigenome = Bio::EnsEMBL::Funcgen::Epigenome->new(
            -dbID               => $eg_id,
            -NAME               => $name,
            -SHORT_NAME         => $short_name,
            -DESCRIPTION        => $desc,
            -GENDER             => $gender,
            -ADAPTOR            => $self,
            -PRODUCTION_NAME    => $production_name,
            -SEARCH_TERMS       => $search_terms_ref,
            -FULL_NAME          => $full_name,
        );

        push @result, $epigenome;

    }
    return \@result;
}

sub _fetch_all_epigenome_ids_having_PeakCalling_by_feature_type_class {

  my $self  = shift;
  my $class = shift;
  
  my $dbc = $self->db->dbc;
  
  use Bio::EnsEMBL::Utils::SqlHelper;
  my $helper = Bio::EnsEMBL::Utils::SqlHelper->new( 
    -DB_CONNECTION => $dbc
  );
  my $epigenome_ids_having_PeakCalling
   = $helper->execute_simple(
    -SQL => '
      select 
        distinct epigenome_id 
      from 
        peak_calling 
        join feature_type using (feature_type_id)
        join epigenome using (epigenome_id)
      where class = ?
     ',
     -PARAMS => [ $class ] 
  );
  return $epigenome_ids_having_PeakCalling
}

sub fetch_all_having_PeakCalling_by_class {

  my $self  = shift;
  my $class = shift;
  
  my $feature_type_ids 
    = $self->_fetch_all_epigenome_ids_having_PeakCalling_by_feature_type_class($class);

  my $feature_types_having_PeakCalling = [
    map {
      $self->fetch_by_dbID($_);
    } 
      @$feature_type_ids
  ];
  
  return $feature_types_having_PeakCalling;

}

=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::Epigenome objects
  Example    : $chan_a->store($c1, $c2, $c3);
  Description: Stores Epigenome objects in the database.
  Returntype : None
  Exceptions : None
  Caller     : General
  Status     : At Risks

=cut

sub store {
  my $self = shift;
  my @args = @_;


  my $sth = $self->prepare("
			INSERT INTO epigenome
			(name, short_name, description, gender, production_name, search_terms, full_name)
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

    my $search_terms;
    if (defined $eg->search_terms){
      $search_terms = stringify($eg->search_terms);
    }

    $sth->bind_param( 1, $eg->name,               SQL_VARCHAR );
    $sth->bind_param( 2, $eg->short_name,         SQL_VARCHAR );
    $sth->bind_param( 3, $eg->description,        SQL_VARCHAR );
    $sth->bind_param( 4, $eg->gender,             SQL_VARCHAR );
    $sth->bind_param( 5, $eg->production_name,    SQL_VARCHAR );
    $sth->bind_param( 6, $search_terms);
    $sth->bind_param( 7, $eg->full_name);
	  
	  $sth->execute();
	  $eg->dbID($self->last_insert_id);
	  $eg->adaptor($self);
	}

  return \@args;
}




1;

