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

  Arg [-STATES]  : Arrayref of states
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

  @{$self->{'states'}} = @$states if $states;

  return $self;
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

   my @state = grep(/$status/, @{$self->get_all_states()});
   my $boolean = scalar(@state);#will be 0 or 1 due to table contraints

   return $boolean;
}



#There is a potential to create an obj from scratch which may already exist in the db
#If we add a state to this (obj has not dbID so will not retrieve stored states) 
# and then try and store it, this will result in adding the state to the previously stored obj.
#The behaviour is silent and could cause problems.

#To resolve this the adaptor implementations must throw if we find a matching object
#We must force the user to generate the obj from the db(use recover) rather than from scratch
#to make them aware of the situation.  This is useful to protect objects where we do not want to overwrite previous data
#e.g. experiment, experimental_chip, channel
#For objects which are routinely resued, we must make sure we always try the db first(not just when recover is set)
#Then warn/throw if there are differing attributes

#This is not possible for set objects, but is not a problem as it will just create another set entry rather than overwriting
#All update/store_states methods should be okay so long as we have a dbID first.




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

   #This could miss states in the DB for storables which have been created and had states added
   #but already exist with states in the DB
   #The way to get around this is to throw if we try and store an object without a dbID which matches 
   #something in the DB.
   #Remove func in adaptors(ec and channel only?) to automatically use prestored objects, throw instead if no dbID and matches.
   #force use of recover to retrieve object from DB and then skip to relevant step based on states.
   #Have states => next method hash in Importer/ArrayDefs?

  

   if($self->is_stored($self->adaptor->db()) && ! $self->{'states'}){
     @{$self->{'states'}} = @{$self->adaptor->fetch_all_states($self)};
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

   if($self->adaptor && $self->is_stored($self->adaptor->db()) && ! $self->{'states'}){
     @{$self->{'states'}} = @{$self->adaptor->fetch_all_states($self)};
   }

   push @{$self->{'states'}}, $status;

   return;
}

sub is_displayable{
  my $self = shift;

  return $self->has_status('DISPLAYABLE');
}

#These DBEntry methods are here to enable xrefs to FeatureType
#Should really on be used for FeatureTypes and SetFeatures

=head2 get_all_Gene_DBEntries

  Example    : my @gene_dbentries = @{ $storable->get_all_Gene_DBEntries };
  Description: Retrieves Ensembl Gene DBEntries (xrefs) for this Storable.  
               This does _not_ include the corresponding translations 
               DBEntries (see get_all_DBLinks).

               This method will attempt to lazy-load DBEntries from a
               database if an adaptor is available and no DBEntries are present
               on the transcript (i.e. they have not already been added or 
               loaded).
  Returntype : Listref of Bio::EnsEMBL::DBEntry objects
  Exceptions : none
  Caller     : general
  Status     : at risk

=cut

#change the to ensembl_core when we implement Gene/Transcript/Protein|Translation links on the same external_db

sub get_all_Gene_DBEntries {
  my $self = shift;;
  return $self->get_all_DBEntries('ensembl_core_Gene');
}

=head2 get_all_Transcript_DBEntries

  Example    : my @transc_dbentries = @{ $set_feature->get_all_Transcript_DBEntries };
  Description: Retrieves ensembl Transcript DBEntries (xrefs) for this Storable.  
               This does _not_ include the corresponding translations 
               DBEntries (see get_all_DBLinks).

               This method will attempt to lazy-load DBEntries from a
               database if an adaptor is available and no DBEntries are present
               on the Storable (i.e. they have not already been added or 
               loaded).
  Returntype : Listref of Bio::EnsEMBL::DBEntry objects
  Exceptions : none
  Caller     : general
  Status     : at risk

=cut

#change the to ensembl_core when we implement Gene/Transcript/Protein|Translation links on the same external_db

sub get_all_Transcript_DBEntries {
  my $self = shift;
  return $self->get_all_DBEntries('ensembl_core_Transcript');
}


=head2 get_all_DBEntries

  Arg[1]     : string - External DB name e.g. ensembl_core_Gene
  Arg[2]     : string - External DB type 
  Example    : my @dbentries = @{ $set_feature->get_all_DBEntries };
  Description: Retrieves DBEntries (xrefs) for this SetFeature.  
               This does _not_ include the corresponding translations 
               DBEntries (see get_all_DBLinks).

               This method will attempt to lazy-load DBEntries from a
               database if an adaptor is available and no DBEntries are present
               on the SetFeature (i.e. they have not already been added or 
               loaded).
  Returntype : Listref of Bio::EnsEMBL::DBEntry objects
  Exceptions : none
  Caller     : general, get_all_DBLinks
  Status     : Stable - at risk move to storable

=cut


#We could add 3rd arg here which would be xref(info_)type e.g. Gene/Transcript etc.
#Move info_type to ox.linkage_type to sit along side linkage_annotated


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

  #This is not using the caching optimally
  #It seems for naive(ex_db_exp,ex_db_type) queries we create a naive cache
  #This means that further more specific queries will make another query and not use the cache


  my @tables = $self->adaptor->_tables;

  if(!defined $self->{$cache_name} && $self->adaptor()) {

	my @tables = $self->adaptor->_tables;
	@tables = split/_/, $tables[0]->[0];#split annotated_feature
	my $object_type = join('', (map ucfirst($_), @tables));#change to AnnotatedFeature
	
    $self->{$cache_name} = 
      $self->adaptor->db->get_DBEntryAdaptor->_fetch_by_object_type($self->dbID(), $object_type, $ex_db_exp, $ex_db_type);
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
