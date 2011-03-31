#
# Ensembl module for Bio::EnsEMBL::Funcgen::Storable
#

=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
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

=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::Storable;


use Bio::EnsEMBL::Registry;
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


  my ($states, $assoc_ftypes) = rearrange(['STATES', 'ASSOCIATED_FEATURE_TYPES'] ,@_);

  if ($self->adaptor() && (! $self->adaptor->isa("Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor"))){
    throw("Adaptor muct be a valid Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor");
  } 

  #will these break using _new_fast
  #THerefore ResultFeature, Probe and ProbeFeature should not be Funcgen::Storables

  @{$self->{'states'}} = @$states if $states;
  $self->associated_feature_types($assoc_ftypes) if(defined $assoc_ftypes);
  

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

#These DBEntry methods are here to enable xrefs to FeatureType, Probe, ProbeSet & ProbeFeature
#They only work as SetFeature has Storable as the second element of @ISA and Bio::EnsEMBL::Feature
#get_all_DBEntries be incorporated into Bio::EnsEMBL::Storable as generic method
#With the other wrapper methods in the Storables of the non-core APIs?
#Or can these be moved to core as a supra core class?
#i.e. Is common to all non-core APIs but not relevant for core?
#Can we bring these together to stop code propogation?


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

#Need to add optional Transcript/Gene param so we can filter
#Filter here or would be better to restrict in sql query ni DBEntryAdaptor?

sub get_all_Gene_DBEntries {
  my $self = shift;


  #We wouldn't need this if we made the xref schema multi species
  #my $species = $self->adaptor->db->species;
  my $species = Bio::EnsEMBL::Registry->get_alias($self->adaptor->db->species);

  if(!$species){
	throw('You must specify a DBAdaptor -species to retrieve DBEntries based on the external_db.db_name');
  }
  
  #safety in case we get Homo sapiens
  ($species = lc($species)) =~ s/ /_/;
  


  return $self->get_all_DBEntries($species.'_core_Gene');
}

=head2 get_all_Transcript_DBEntries

  Arg[0]     : optional - Bio::EnsEMBL::Transcript to filter DBEntries on.
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
  my ($self, $transcript) = @_;

  #We wouldn't need this if we made the xref schema multi species
  my $species = Bio::EnsEMBL::Registry->get_alias($self->adaptor->db->species);
  #Need to make sure this is latin name
  

  if(!$species){
	throw('You must specify a DBAdaptor -species to retrieve DBEntries based on the external_db.db_name');
  }
  
  #safety in case we get Homo sapiens
  #($species = lc($species)) =~ s/ /_/;
  
  my $dbes = $self->get_all_DBEntries($species.'_core_Transcript');

  #This needs to be moved to the DBEntryAdaptor and restrict the query using the
  #dbprimary_acc

  if($transcript){
	my $sid = $transcript->stable_id;

	#Test for sid here?

	if(ref($transcript) && $transcript->isa('Bio::EnsEMBL::Transcript')){
	  my @dbes;

	  foreach my $dbe(@$dbes){
		if($dbe->primary_id eq $sid){
		  push @dbes, $dbe;
		}
	  }
	  $dbes = \@dbes;
	}
	else{
	  throw('Transcript argument must be a valid Bio::EnsEMBL::Transcript');
	}
  }


  return $dbes;
}


=head2 get_all_UnmappedObjects

  Example    : my @uos = @{$storable->get_all_UnmappedObjects };
  Description: Retrieves UnamappedObjects for this Storable.
  Returntype : arrayref of Bio::EnsEMBL::UnmappedObject objects
  Exceptions : none
  Caller     : general
  Status     : At risk - move to Bio::Ensembl::Storable?

=cut

sub get_all_UnmappedObjects{
  my $self = shift;
  #Do we want to add external_db or analysis param here?

  my $class = ref($self);
  $class =~ s/.*:://;

  return $self->adaptor->db->get_UnmappedObjectAdaptor->fetch_all_by_object_type_id($class, $self->dbID);
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


  if( (! defined $self->{$cache_name}) && $self->adaptor() ){
  
	my @tables = $self->adaptor->_tables;
	@tables = split/_/, $tables[0]->[0];#split annotated_feature
	my $object_type = join('', (map ucfirst($_), @tables));#change to AnnotatedFeature
	
    $self->{$cache_name} = 
      $self->adaptor->db->get_DBEntryAdaptor->_fetch_by_object_type($self->dbID(), $object_type, $ex_db_exp, $ex_db_type);
  }
  elsif( ! defined $self->{$cache_name} ){
	throw('You must have set and adaptor to be able to get_all_DBEntries');
  }


  $self->{$cache_name} ||= [];

  return $self->{$cache_name};
}


=head2 add_DBEntry

  Arg [1]    : Bio::EnsEMBL::DBEntry $dbe
               The dbEntry to be added
  Example    : my $dbe = Bio::EnsEMBL::DBEntry->new(...);
               $transcript->add_DBEntry($dbe);
  Description: Associates a DBEntry with this object. Note that adding
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


=head2 associated_feature_types

  Example    : my @associated_ftypes = @{$feature->associated_feature_types()};
  Description: Getter/Setter for other associated FeatureTypes.
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen:FeatureType objects
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub associated_feature_types{
  my ($self, $ftypes) = @_;
  
  #Lazy load as we don't want to have to do a join on all features when most will not have any

 
  if(defined $ftypes){

	if(ref($ftypes) eq 'ARRAY'){

	  foreach my $ftype(@$ftypes){
	
		if( ! $ftype->isa('Bio::EnsEMBL::Funcgen::FeatureType') ){
		  throw('You must pass and ARRAYREF of stored Bio::EnsEMBL::Funcgen::FeatureType objects');
		}
		#test is stored in adaptor
	  }

	  if(defined $self->{'associated_feature_types'}){
		warn('You are overwriting associated feature types');
		#we could simply add the new ones and make them NR.
	  }

	  $self->{'associated_feature_types'} = $ftypes;
	}
	else{
	  throw('You must pass and ARRAYREF of stored Bio::EnsEMBL::Funcgen::FeatureType objects');
	}
  }


  if(! defined $self->{'associated_feature_types'}){
	#This will fail if we have not stored yet

	if(defined $self->adaptor){
	  $self->{'associated_feature_types'} = $self->adaptor->db->get_FeatureTypeAdaptor->fetch_all_by_association($self);
	}

  }


  #This has the potential to return undef, or an arrayref which may be empty.
  return $self->{'associated_feature_types'};
}







1;
