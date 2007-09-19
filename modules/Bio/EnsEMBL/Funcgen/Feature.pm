# Ensembl module for Bio::EnsEMBL::Funcgen::Feature
#
# Copyright (c) 2003 Ensembl
#


=head1 NAME

Bio::EnsEMBL::Funcgen::Feature - Ensembl specific sequence feature.

=head1 SYNOPSIS

    my $feat = new Bio::EnsEMBL::Feature(-start   => 100,
                                       -end     => 220,
                                       -strand  => -1,
                                       -slice   => $slice
                                       -analysis => $analysis
                                      );

    my $start  = $feat->start;
    my $end    = $feat->end;
    my $strand = $feat->strand;

    #move the feature to the chromosomal coordinate system
    $feature = $feature->transform('chromosome');

    #move the feature to a different slice (possibly on another coord system)
    $feature = $feature->transfer($new_slice);

    #project the feature onto another coordinate system possibly across
    #boundaries:
    @projection = @{$feature->project('contig')};

    #change the start, end, and strand of the feature in place
    $feature->move($new_start, $new_end, $new_strand);

=head1 DESCRIPTION

This is a simple wrapper method for the core Feature class to contain generic
Funcgen Feature DBEntry methods.

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=head1 METHODS

=cut


use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::Feature;

use Bio::EnsEMBL::Feature;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning);
use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Feature);



=head2 get_all_DBEntries

  Example    : my @dbentries = @{ $gene->get_all_DBEntries };
  Description: Retrieves DBEntries (xrefs) for this transcript.  
               This does _not_ include the corresponding translations 
               DBEntries (see get_all_DBLinks).

               This method will attempt to lazy-load DBEntries from a
               database if an adaptor is available and no DBEntries are present
               on the transcript (i.e. they have not already been added or 
               loaded).
  Returntype : Listref of Bio::EnsEMBL::DBEntry objects
  Exceptions : none
  Caller     : get_all_DBLinks, TranscriptAdaptor::store
  Status     : Stable

=cut

sub get_all_DBEntries {
  my $self = shift;
  my $ex_db_exp = shift;
  my $ex_db_type = shift;

  my $cache_name = "dbentries";

  if(defined($ex_db_exp)){
    $cache_name .= $ex_db_exp;
  }
  if(defined($ex_db_type)){
    $cache_name .= $ex_db_type;
  }

  #Need to add tests for valid objects for xrefs

  # if not cached, retrieve all of the xrefs for this gene

  my @tables = $self->adaptor->_tables;

  if(!defined $self->{$cache_name} && $self->adaptor()) {
    $self->{$cache_name} = 
      $self->adaptor->db->get_DBEntryAdaptor->_fetch_by_object_type($self->dbID(), $tables[0]->[0], $ex_db_exp, $ex_db_type);
  }

  $self->{$cache_name} ||= [];

  return $self->{$cache_name};
}


=head2 add_DBEntry

  Arg [1]    : Bio::EnsEMBL::DBEntry $dbe
               The dbEntry to be added
  Example    : my $dbe = Bio::EnsEMBL::DBEntery->new(...);
               $transcript->add_DBEntry($dbe);
  Description: Associates a DBEntry with this transcript. Note that adding
               DBEntries will prevent future lazy-loading of DBEntries for this
               storable (see get_all_DBEntries).
  Returntype : none
  Exceptions : thrown on incorrect argument type
  Caller     : general
  Status     : Stable

=cut

sub add_DBEntry {
  my $self = shift;
  my $dbe = shift;

  unless($dbe && ref($dbe) && $dbe->isa('Bio::EnsEMBL::DBEntry')) {
    throw('Expected DBEntry argument');
  }

  $self->{'dbentries'} ||= [];
  push @{$self->{'dbentries'}}, $dbe;
}



1;
