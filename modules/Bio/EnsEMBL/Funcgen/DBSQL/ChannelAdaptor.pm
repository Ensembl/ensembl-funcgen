#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::ChannelAdaptor
#
# You may distribute this module under the same terms as Perl itself

=head1 NAME

Bio::EnsEMBL::Funcgen::DBSQL::ChannelAdaptor - A database adaptor for fetching and
storing Funcgen Channel objects.

=head1 SYNOPSIS

my $chan_a = $db->get_ChannelAdaptor();

my @channels = @{$chan_a->fetch_all_by_ExperimentalChip($exp)};


=head1 DESCRIPTION

The ChannelAdaptor is a database adaptor for storing and retrieving
Funcgen Channel objects.

=head1 AUTHOR

This module was created by Nathan Johnson.

This module is part of the Ensembl project: http://www.ensembl.org/

=head1 CONTACT

Post comments or questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::DBSQL::ChannelAdaptor;

use Bio::EnsEMBL::Utils::Exception qw( warning throw );
use Bio::EnsEMBL::Funcgen::Channel;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;

use vars qw(@ISA);


#May need to our this?
@ISA = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

=head2 fetch_by_type_experimental_chip_id

  Arg [1]    : string - type
  Arg [1]    : int - dbID of Experiment
  Example    : my $chan = $chan_a->fetch_by_type_experimental_chip_dbID('EXPERIMENTAL', $ec_dbid);
  Description: Does what it says on the tin
  Returntype : Bio::EnsEMBL::Funcgen::Channel object
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub fetch_by_type_experimental_chip_dbID {
    my ($self, $type, $ec_dbid) = @_;

	throw("Must specify a valid experiemntal dbID($ec_dbid) and type($type)") if(! $ec_dbid || ! $type);


	my $sth = $self->prepare("
		SELECT c.channel_id
		FROM channel c
		WHERE c.experimental_chip_id = ?
        AND c.type = ?
	");

	$sth->bind_param(1, $ec_dbid,   SQL_INTEGER);
	$sth->bind_param(2, $type,      SQL_VARCHAR);

	#can we do a generic fetch here?


	$sth->execute();
	my ($chan_id) = $sth->fetchrow();


	return $self->fetch_by_dbID($chan_id) if $chan_id;
}


=head2 fetch_by_dye_experimental_chip_id

  Arg [1]    : string - dye
  Arg [1]    : int - dbID of Experiment
  Example    : my $chan = $chan_a->fetch_by_type_experimental_chip_dbID('Cy5', $ec_dbid);
  Description: Does what it says on the tin
  Returntype : Bio::EnsEMBL::Funcgen::Channel object
  Exceptions : Throws is experiment dbID or dye not passed
  Caller     : General
  Status     : At Risk

=cut

sub fetch_by_dye_experimental_chip_dbID {
    my $self = shift;
    my $dye = shift;
    my $ec_dbid = shift;

    throw("not yet impemented");

    my ($chan_id, @results);

    throw("Must specify an experiemntal dbID") if(! $ec_dbid);
    #Need to validate dye against VendorDefs?  Or leave and just return whatever is in DB e.g. nothing if dye name is wrong.


    my $sth = $self->prepare("
		SELECT c.channel_id
		FROM experimental_chip ec, channel c
		WHERE c.experimental_chip_id = ec.experimental_chip_id
        AND ec.experimental_chip_id = $ec_dbid
	");



	#can we do a generic fetch here?


	$sth->execute();


	while ($chan_id = $sth->fetchrow()){
		push @results, $self->fetch_by_dbID($chan_id);
	}

	return \@results;
}


=head2 fetch_all_by_experimental_chip_dbID

  Arg [1]    : int - dbID of Experiment
  Example    : my @chans = @{$ec_a->fetch_all_by_experimental_chip_dbID($ac_dbid);
  Description: Does what it says on the tin
  Returntype : Listref of Bio::EnsEMBL::Funcgen::Channel object
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_experimental_chip_dbID {
    my $self = shift;
    my $ec_dbid = shift;

	my ($chan_id, @results);

	throw("Must specify an experiemntal dbID") if(! $ec_dbid);


	my $sth = $self->prepare("
		SELECT c.channel_id
		FROM experimental_chip ec, channel c
		WHERE c.experimental_chip_id = ec.experimental_chip_id
        AND ec.experimental_chip_id = $ec_dbid
	");


	#can we do a generic fetch here?


	$sth->execute();


	while ($chan_id = $sth->fetchrow()){
		push @results, $self->fetch_by_dbID($chan_id);
	}

	return \@results;
}

=head2 fetch_all_by_ExperimentalChip

  Arg [1]    : Bio::EnsEMBL::Funcgen::ExperimentalChip
  Example    : my @chans = @{$ec_a->fetch_all_by_ExperimentalChip($echip);
  Description: Returns all channels associated with a given ExperimentalChip
  Returntype : Listref of Bio::EnsEMBL::Funcgen::Channel objects
  Exceptions : Throws if no ExperimentalChip defined
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_ExperimentalChip{
  my ($self, $exp) = @_;
  throw("Must provide an ExperimentChip object") if(! $exp->isa('Bio::EnsEMBL::Funcgen::ExperimentalChip'));
  return $self->fetch_all_by_experimental_chip_dbID($exp->dbID());
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
	
	return ['channel', 'c'];
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
	
	return qw( c.channel_id c.experimental_chip_id c.sample_id c.cell_line_id c.dye c.type c.description);
}

=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates Channel objects from an executed DBI statement
			   handle.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::Channel objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
	my ($self, $sth) = @_;
	
	my (@result, $ec_id, $chan_id, $sample_id, $cell_line_id, $type, $dye, $desc);
	
	$sth->bind_columns(\$chan_id, \$ec_id, \$sample_id, \$cell_line_id, \$dye, \$type, \$desc);
	
	while ( $sth->fetch() ) {
	  my $chan = Bio::EnsEMBL::Funcgen::Channel->new(
							 -dbID                 => $chan_id,
							 -EXPERIMENTAL_CHIP_ID => $ec_id,
							 -SAMPLE_ID            => $sample_id,
							 -CELL_LINE_ID         => $cell_line_id,
							 -TYPE                 => $type,
							 -DYE                  => $dye,
							 -DESCRIPTION          => $desc,
							 -ADAPTOR              => $self,
							);
	  
	  push @result, $chan;
	  
	}
	return \@result;
}



=head2 store

  Args       : List of Bio::EnsEMBL::Funcgen::Channel objects
  Example    : $chan_a->store($c1, $c2, $c3);
  Description: Stores given Channel objects in the database. Should only be
               called once per array because no checks are made for duplicates.
			   Sets dbID and adaptor on the objects that it stores.
  Returntype : None
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub store {
  my $self = shift;
  my @args = @_;
  
  my ($sarray);
  
  my $sth = $self->prepare("
			INSERT INTO channel
			(experimental_chip_id, sample_id, cell_line_id, dye, type, description)
			VALUES (?, ?, ?, ?, ?, ?)");
    
  
  
  foreach my $chan (@args) {
    if ( ! $chan->isa('Bio::EnsEMBL::Funcgen::Channel') ) {
      warning('Can only store Channel objects, skipping $chan');
      next;
    }
    
    if (!( $chan->dbID() && $chan->adaptor() == $self )){
      
      
      my $s_chan = $self->fetch_by_type_experimental_chip_dbID($chan->type(), $chan->experimental_chip_id());
	
	
	if(! $s_chan){
	  $sth->bind_param(1, $chan->experimental_chip_id(),  SQL_INTEGER);
	  $sth->bind_param(2, $chan->sample_id(),             SQL_VARCHAR);
	  $sth->bind_param(3, $chan->cell_line_id(),          SQL_INTEGER);
	  $sth->bind_param(4, $chan->dye() ,                  SQL_VARCHAR);
	  $sth->bind_param(5, $chan->type(),                  SQL_VARCHAR);
	  $sth->bind_param(6, $chan->description(),           SQL_VARCHAR);
	  
	  $sth->execute();
	  my $dbID = $sth->{'mysql_insertid'};
	  $chan->dbID($dbID);
	  $chan->adaptor($self);
	}
	else{
	  #do some status checks here, check IMPORTED
	  #Need to account for recover in Importer?
	  $chan = $s_chan;

	  my @states = @{$self->db->fetch_all_states('channel', $chan->dbID())};

	  #need better id than dbID?
	  warn("Using previously stored Channel (".$chan->experimental_chip_id().":".$chan->type().") with states\t@states\n"); 
	}
      }
    }

    return \@args;

}


=head2 list_dbIDs

  Args       : None
  Example    : my @array_ids = @{$ec_a->list_dbIDs()};
  Description: Gets an array of internal IDs for all ExperimentalChip objects in the
               current database.
  Returntype : List of ints
  Exceptions : None
  Caller     : ?
  Status     : At risk

=cut

sub list_dbIDs {
    my ($self) = @_;
	
    return $self->_list_dbIDs('experimental_chip');
}



1;

