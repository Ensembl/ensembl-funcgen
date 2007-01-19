#
# Ensembl module for Bio::EnsEMBL::Funcgen::Storable
#

=head1 NAME

Bio::EnsEMBL::Funcgen::Storable

=head1 SYNOPSIS

  my $dbID = $storable_object->dbID();
  my $adaptor = $storable_object->adaptor();
  if($storable_object->is_stored($db_adaptor))) {
    
  }
=head1 DESCRIPTION

This is a simple wrapper class to provide convenience methods for the StorableAdaptor.
Only get type methods have been implemented here to avoid obfuscating DB writes which 
should only be done by the specific 'Storable'Adaptors.

=head1 AUTHOR - Nathan Johnson

=head1 CONTACT

Post questions to the Ensembl development list B<ensembl-dev@ebi.ac.uk>

=head1 METHODS

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::Storable;



use Bio::EnsEMBL::Utils::Exception qw(throw warning);
use Bio::EnsEMBL::Utils::Argument qw(rearrange);
use Bio::EnsEMBL::Storable;

use vars qw(@ISA);

@ISA = qw(Bio::EnsEMBL::Storable);

=head2 new

  Arg [-ADAPTOR] : Bio::EnsEMBL::DBSQL::BaseAdaptor
  Arg [-dbID]    : database internal id
  Example        : none 
  Caller         : internal calls
  Description    : create a new Storable object 
  Returntype     : Bio::EnsEMBL::Storable
  Exceptions     : Adaptor not a Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor
  Status         : Stable

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;
	
  my $self = $class->SUPER::new(@_);


  my ($states) = rearrange(['STATES'],@_);

  if ($self->adaptor() && (! $self->adaptor->isa("Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor"))){
    throw("Adaptor muct be a valid Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor");
  } 

  #will these break using _new_fast
  #THerefore ResultFeature, Probe and ProbeFeature should not be Funcgen::Storables


  $self->{'states'} = ();
  $self->{'states'} = @$states if $states;


  return $self;
}


=head2 _status_adaptor

  Example    : if($self->_status_adaptor->has_state($self, $status){ ... }
  Description: Internal convenience accessor to StatusAdaptor
  Returntype : Bio::EnsEMBL::Funcgen::DBSQL::StatusAdaptor
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::Storable
  Status     : At risk

=cut



sub _status_adaptor{
   my ($self) = @_;

   $self->{'_status_adaptor'} ||= $self->adaptor->db->get_StatusAdaptor();
   
   return $self->{'_status_adaptor'};
}


=head2 has_status

  Arg [1]    : string - status e.g. IMPORTED, DISPLAYABLE
  Example    : if($experimental_chip->has_status('IMPORTED'){ ... skip import ... };
  Description: Tests whether storable has a given status
  Returntype : BOOLEAN
  Exceptions : Throws if not status is provided
  Caller     : general
  Status     : At risk

=cut



sub has_status{
   my ($self, $status) = @_;

   throw("Must provide a status to check") if ! $status;

   my @state = grep(/$status/, @{$self->{'states'}});
   my $boolean = scalar(@state);#will be 0 or 1 due to table contraints

   return $boolean;
}



=head2 get_all_states

  Example    : my @ec_states = @{$experimental_chip->get_all_states()};
  Description: Retrieves all states from DB and merges with current states array
  Returntype : LISTREF
  Exceptions : None
  Caller     : general
  Status     : At risk

=cut



sub get_all_states{
   my ($self) = @_;

   my %states;

   #overwriting the cache is prevented by add_states getting states from db if already stored

   if($self->is_stored() && ! $self->{'states'}){
     $self->{'states'} = @{$self->_status_adaptor->fetch_all_states($self)};
   }

   return $self->{'states'};
}


=head2 add_status

  Example    : $ec->add_state('DISPLAYABLE');
  Description: Adds a state to a new or previously stored Storable
  Returntype : None
  Exceptions : Throws if no status supplied
  Caller     : general
  Status     : At risk

=cut



sub add_status{
   my ($self, $status) = @_;

   throw("Must pass a status to add e.g. 'DISPLAYABLE'") if ! $status;



   #this does not resolve the problem!!???
   #can add a status to an unstored object which 

   if($self->is_stored() && ! $self->{'states'}){
     $self->{'states'} = @{$self->adaptor->fetch_all_states($self)};
   }

   push @{$self->{'states'}}, $status;

   return;
}

sub is_displayable{
  my $self = shift;

  return $self->has_status('DISPLAYABLE');
}

#Need to update adaptor store methods for ExperimentalChip, FeatureSet, ResultSet and Channel?
#Add store_states, do we have a Funcgen::BaseFeatureAdaptor


1;
