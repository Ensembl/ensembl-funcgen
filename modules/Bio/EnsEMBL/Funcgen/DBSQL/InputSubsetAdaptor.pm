#
# Ensembl module for Bio::EnsEMBL::Funcgen::
#

=head1 LICENSE

  Copyright (c) 1999-2013 The European Bioinformatics Institute and
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

Bio::EnsEMBL::Funcgen::DBSQL::InputSubsetAdaptor - A module to


=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::DBSQL::InputSubsetAdaptor
my $iss_a = $db->get_InputSubsetAdaptor;


=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Funcgen::DBSQL::InputSubsetAdaptor;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(dump_data);
use Bio::EnsEMBL::Funcgen::InputSubset;
use Bio::EnsEMBL::Funcgen::DBSQL::SetAdaptor;

use base qw(Bio::EnsEMBL::Funcgen::DBSQL::SetAdaptor);



=head2 fetch_all_by_InputSet

  Arg [1]    : Bio::EnsEMBL::Funcgen::InputSet
  Example    : my @i_subsets = @{$is_adaptopr->fetch_all_by_InputSet($type)};
  Description: Retrieves InputSubset objects from the database based on InputSet
  Returntype : Listref of Bio::EnsEMBL::Funcgen::InputSubset objects
  Exceptions : Throws if arg is not a valid InputSet
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_InputSet {
  my ($self, $iset) = @_;

  my $params = {constraints => {input_sets => [$iset]}};

#  $self->bind_param_generic_fetch($ftype->dbID,  SQL_INTEGER);
#  my $objs = $self->generic_fetch($constraint);
#  $self->reset_true_tables;

	my $objs = $self->generic_fetch($self->compose_constraint_query($params));

  $self->reset_true_tables;

  return($objs);
}

=head2 fetch_all_by_Experiment

  Arg [1]    : Bio::EnsEMBL::Funcgen::Experiment
  Example    : my @i_subsets = @{$ex_adaptopr->fetch_all_by_Experiment($type)};
  Description: Retrieves InputSubset objects from the database based on Experiment
  Returntype : Listref of Bio::EnsEMBL::Funcgen::InputSubset objects
  Exceptions : Throws if arg is not a valid Experiment
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_Experiments {
  my ($self, $exps) = @_;

  my $params = {constraints => {experiments => [$exps]}};
	return $self->generic_fetch($self->compose_constraint_query($params));

}

=head2 fetch_by_name

  Arg [1]    : string - InputSubset name
  Example    : my $iss = $iss_a->fetch_by_name('Iss_name');
  Description: Retrieves a InputSubset object which matches the passed name
  Returntype : Bio::EnsEMBL::Funcgen::InputSubset
  Exceptions : Throws if no name defined or if more than one returned
  Caller     : General
  Status     : At risk

=cut

sub fetch_all_by_name {
  my ($self, $name) = @_;

  throw("Need to specify a name argument") if (! defined $name);
  $self->bind_param_generic_fetch($name, SQL_VARCHAR);
  return $self->generic_fetch("name = ?")->[0];

}

=head2 fetch_by_name_and_experiment

  Arg [1]    : string - InputSubset name
  Arg [2]    : Bio::EnsEMBL::Funcgen::Experiment
  Example    : my $iss = $iss_a->fetch_by_name_and_experiment('Iss_name', $exp);
  Description: Retrieves a InputSubset object which matches the passed name and experiment
  Returntype : Bio::EnsEMBL::Funcgen::InputSubset
  Exceptions : Throws if no name defined or if more than one returned
  Caller     : General
  Status     : At risk

=cut

sub fetch_by_name_and_experiment {
  my ($self, $name, $exp) = @_;

  my $params = {
    constraints => {
      experiments => [$exp],
      name        => $name
    }};
	return $self->generic_fetch($self->compose_constraint_query($params));
}


=head2 _columns

  Args       : None
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Returns a list of columns to use for creating the object
               and queries.
  Returntype : List of strings
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _columns {
  my ($self) = @_;

  return qw(
      iss.input_subset_id
      iss.experiment_id
      iss.archive_id
      iss.display_url
      iss.is_control
      iss.name
      iss.replicate
      );
}



=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates Array objects from an executed DBI statement
               handle.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::InputSubset objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
  my ($self, $sth) = @_;

  my (@result,
      $iss_id,
      $exp_id,
      $archive_id,
      $display_url,
      $is_control,
      $name,
      $replicate,
     );

#  my $is_adaptor = $self->db->get_InputSetAdaptor;

  #Fetching the InputSet will need to be removed once we implement the lin ktable as there maybe >!
  #so we would need a hash of input_sets => {dbDB1 => $iset1, dbID2 => $iset2}
  #we're initially the values are undef, and lazy loaded on request.

  $sth->bind_columns(
      \$iss_id,
      \$exp_id,
      \$archive_id,
      \$display_url,
      \$is_control,
      \$name,
      \$replicate,
      );

  my $exp_adaptor = $self->db->get_ExperimentAdaptor;

  while($sth->fetch()){
    my $exp = $exp_adaptor->fetch_by_dbID($exp_id);
    if(! defined $exp){
      throw("Could not fetch linked experiment (dbID: $exp_id) for InputSubset (dbID: $iss_id) ");
    }

    push @result, Bio::EnsEMBL::Funcgen::InputSubset->new (
       -dbID         => $iss_id,
       -EXPERIMENT   => $exp,
       -ARCHIVE_ID   => $archive_id,
       -DISPLAY_URL  => $display_url,
       -IS_CONTROL   => $is_control,
       -NAME         => $name,
       -REPLICATE    => $replicate,
      );
  }
  return(\@result);
}

=head2 store

  Args       : Bio::EnsEMBL::Funcgen::InputSet
  Example    : $iss_adaptor->store(\@e_subsets);
  Description: Stores or updates previously stored InputSubset objects in the database.
  Returntype : None
  Exceptions : Throws if a List of InputSet objects is not provided or if
               an analysis is not attached to any of the objects
  Caller     : General
  Status     : At Risk

=cut


sub store{
  my ($self, @subsets) = @_;

  if(scalar(@subsets == 0)){
    throw("Must provide a list of InputSubset objects")
  }

  my $db = $self->db();

  my $sth = $self->prepare("
        INSERT INTO
          input_subset (
            experiment_id,
            archive_id,
            display_url,
            is_control,
            name,
            replicate
        )
        VALUES (?, ?, ?, ?, ? ,?)
    ");


  #Store and set all previously unstored table_ids
  foreach my $subset(@subsets){
    if( ! ref $subset || ! $subset->isa('Bio::EnsEMBL::Funcgen::InputSubset') ) {
      throw("InputSubset object required, not '".ref($subset) ."'");
    }

    if($subset->is_stored($db)){
      throw( "InputSubset '" . $subset->name . "' already stored in the DB.".
            "InputSubsetAdaptor does not support updating yet.");
    }

    $sth->bind_param(1, $subset->experiment->dbID, SQL_INTEGER);
    $sth->bind_param(2, $subset->archive_id,       SQL_VARCHAR);
    $sth->bind_param(3, $subset->display_url,      SQL_VARCHAR);
    $sth->bind_param(4, $subset->is_control,       SQL_INTEGER);
    $sth->bind_param(5, $subset->name,             SQL_VARCHAR);
    $sth->bind_param(6, $subset->replicate,        SQL_INTEGER);
    $sth->execute();
# Core utils fetch last ID, platform independant
    $subset->dbID($self->last_insert_id);
    $subset->adaptor($self);
    $self->store_states($subset);
  }

  return \@subsets;
}


sub _final_clause {
  return ( 'group by iss.input_subset_id');
}

sub _constrain_experiments {
  my ($self, $exps) = @_;

  my $constraint = ' iss.experiment_id IN ('.
    join(', ', @{$self->db->are_stored_and_valid('Bio::EnsEMBL::Funcgen::Experiment', $exps, 'dbID')}).')';

  #{} = not futher constraint conf
  return ($constraint, {});

}

sub _constrain_input_sets {
  my ($self, $isets) = @_;

  my $constraint_conf = {
    tables => [
      ['input_set_input_subset', 'isiss']
      ]};

  my $constraint = '
    iss.input_subset_id = isiss.input_subset_id AND
    isiss.input_set_id  IN ('.
    join(', ', @{$self->db->are_stored_and_valid('Bio::EnsEMBL::Funcgen::InputSet', $isets, 'dbID')}).')';

#  return ($constraint);
return ($constraint, $constraint_conf);
}

sub _constrain_archive_ids {
  my ($self, $archive_ids) = @_;

  if ( (ref($archive_ids) ne 'ARRAY') || scalar(@$archive_ids) == 0 ) {
    throw('Must pass an arrayref of archive IDs');
  }
  my $constraint;

  #{} = not futher constraint conf
  return (' iss.archive_id IN '.join(', ', (map {uc($_)} @$archive_ids)) , {});
}

sub _true_tables {
  return ([ 'input_subset',  'iss' ]);
}

1;

