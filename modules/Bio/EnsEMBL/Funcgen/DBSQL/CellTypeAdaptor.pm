#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::CellTypeAdaptor
#

=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.


=head1 NAME

Bio::EnsEMBL::Funcgen::DBSQL::CellTypeAdaptor - A database adaptor for fetching and
storing Funcgen CellType objects.

=head1 SYNOPSIS

my $ct_adaptor = $efgdba->get_CellTypeAdaptor();

my $cell_type = $ct_adaptor->fetch_by_name("HeLa-S3");


=head1 DESCRIPTION

The CellTypeAdaptor is a database adaptor for storing and retrieving
Funcgen CellType objects.

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::DBSQL::CellTypeAdaptor;

use Bio::EnsEMBL::Utils::Exception qw( warning throw );
use Bio::EnsEMBL::Funcgen::CellType;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;

use vars qw(@ISA);


#May need to our this?
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

=head2 fetch_by_name

  Arg [1]    : string - name of CellType
  Arg [1]    : optional string - class of CellType
  Example    : my $ct = $ct_adaptor->fetch_by_name('HeLa');
  Description: Retrieves CellType objects by name.
  Returntype : Bio::EnsEMBL::Funcgen::CellType object
  Exceptions : Throws no name given
  Caller     : General
  Status     : At risk

=cut

sub fetch_by_name{
  my ($self, $name) = @_;

  throw("Must specify a CellType name") if(! $name);

  my $constraint = "ct.name ='$name'";

  my @ctype = @{$self->generic_fetch($constraint)};
  #name is unique so we should only have one

  return $ctype[0];
}


#fetch_all_by_efo_id
#fetch_all_by_tissue
#fetch_all_by_lineage

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
	  ['cell_type', 'ct'],
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
	
  return qw( ct.cell_type_id ct.name ct.display_label ct.description ct.gender ct.efo_id ct.tissue);
  #type/class = enum cell, cell line, tissue
}

=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates Channel objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::CellType objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
	my ($self, $sth) = @_;
	
	my (@result, $ct_id, $name, $dlabel, $desc, $gender, $efo_id, $tissue);
	
	$sth->bind_columns(\$ct_id, \$name, \$dlabel, \$desc, \$gender, \$efo_id, \$tissue);
	
	while ( $sth->fetch() ) {
		my $ctype = Bio::EnsEMBL::Funcgen::CellType->new(
														 -dbID          => $ct_id,
														 -NAME          => $name,
														 -DISPLAY_LABEL => $dlabel,
														 -DESCRIPTION   => $desc,
														 -GENDER        => $gender,
														 -EFO_ID        => $efo_id,
														 -ADAPTOR       => $self,
														);
	  
		push @result, $ctype;
	  
	}
	return \@result;
}



=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::CellType objects
  Example    : $chan_a->store($c1, $c2, $c3);
  Description: Stores CellType objects in the database.
  Returntype : None
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub store {
  my $self = shift;
  my @args = @_;
  
  
  my $sth = $self->prepare("
			INSERT INTO cell_type
			(name, display_label, description, gender, efo_id, tissue)
			VALUES (?, ?, ?, ?, ?, ?)");
    
  
  
  foreach my $ct (@args) {
	  if ( ! $ct->isa('Bio::EnsEMBL::Funcgen::CellType') ) {
		  warning('Can only store CellType objects, skipping $ct');
		  next;
	  }
	  
	  if ( $ct->dbID() && $ct->adaptor() == $self ){
		  warn("Skipping previously stored CellType dbID:".$ct->dbID().")");
		  next;
	  }
	  
	  
	  $sth->bind_param(1, $ct->name,           SQL_VARCHAR);
	  $sth->bind_param(2, $ct->display_label,  SQL_VARCHAR);
	  $sth->bind_param(3, $ct->description,    SQL_VARCHAR);
	  $sth->bind_param(4, $ct->gender,         SQL_VARCHAR);
	  $sth->bind_param(5, $ct->efo_id,         SQL_VARCHAR);
	  $sth->bind_param(6, $ct->tissue,         SQL_VARCHAR);
	  $sth->execute();
	  $ct->dbID($sth->{'mysql_insertid'});
	  $ct->adaptor($self);
	}

  return \@args;
}




1;

