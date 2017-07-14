#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::ProbeAdaptor
#

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

Bio::EnsEMBL::Funcgen::DBSQL::ProbeAdaptor - A database adaptor for fetching and
storing Probe objects.

=head1 SYNOPSIS

my $probe_adaptor = $db->get_ProbeAdaptor();

my $probe = $probe_adaptor->fetch_by_array_probe_probeset('Array-1', 'Probe-1');

=head1 DESCRIPTION

The ProbeAdaptor is a database adaptor for storing and retrieving
Probe objects.

=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::Probe
Bio::EnsEMBL::Funcgen::ProbeFeature
Bio::EnsEMBL::Funcgen::ProbeSet
Bio::EnsEMBL::Funcgen::ArrayChip
Bio::EnsEMBL::Funcgen::Array

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::ProbeAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Utils::Exception qw( deprecate );
use Bio::EnsEMBL::Funcgen::Probe;
use Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor;#DBI sql_types import
use feature qw(say);
use base qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor);

use DBI;

# to standard output
# DBI->trace(1);

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  
  $self->{_array_cache}     = {};
  $self->{_probe_set_cache} = {};
  
  return $self;
}

=head2 fetch_by_array_probe_probeset_name

  Arg [1]    : string - name of array
  Arg [2]    : string - name of probe
  Arg [3]    : (optional) string - name of probeset
  Example    : my $probe = $opa->fetch_by_array_probeset_probe('Array-1', 'Probe-1', 'ProbeSet-1');
  Description: Returns a probe given a combination of array name, probeset and
               probe name. This will uniquely define an Affy probe. Only one
			   probe is ever returned.
  Returntype : Bio::EnsEMBL::Funcgen::Probe
  Exceptions : throws if array or probe name not defined
  Caller     : General
  Status     : At Risk - rename to fetch_by_probe_array_probeset_name?

=cut

#This does not currently capture on plate replicate probes with different names
#Only returns the record corresponding to the given name and not the other replicate name
################# ----------------------- Continue
sub fetch_by_array_probe_probeset_name {
	my ($self, $array_name, $probe_name, $probeset_name) = @_;

	if(! (defined $array_name && defined $probe_name)){
	  throw('You must provide at least and array and probe name');
	}

	my $tables = 'probe p,  array a';
	$tables .= ', probe_set ps' if defined $probeset_name;

	my $sql = "SELECT distinct(p.probe_id) FROM $tables WHERE a.name=? and a.array_id=p.array_id and p.name=?";
	$sql .= ' AND p.probe_set_id=ps.probe_set_id and ps.name=?' if defined $probeset_name;

	my $sth = $self->db->dbc->prepare($sql);
	$sth->bind_param(1, $array_name,    SQL_VARCHAR);
	$sth->bind_param(2, $probe_name,    SQL_VARCHAR);
	$sth->bind_param(3, $probeset_name, SQL_VARCHAR) if defined $probeset_name;
	$sth->execute;

	#This should only return one result
	#The only possible way this would not return one result
	#is if an identically named array(:probeset):probe which had a different sequence
	#As Import array would separate these based on the sequence hash
	my ($dbid) = $sth->fetchrow_array;

	return (defined $dbid) ? $self->fetch_by_dbID($dbid) : undef;
}


=head2 fetch_all_by_external_name

  Arg [1]    : Ensembl stable Transcript ID
  Example    : $probeAdaptor->fetch_all_by_external_name
  Description: 
               
  Returntype : 
  Caller     : General
  Status     : Deprecated

=cut

sub fetch_all_by_external_name {
  my $self = shift;
  my $transcript_stable_id = shift;
  deprecate( "fetch_all_by_external_namehas been deprecated and will be removed in Ensembl release 92.");
  return $self->fetch_all_by_transcript_stable_id($transcript_stable_id);
}

=head2 fetch_all_by_sequence

  Arg [1]    : string sequence
  Example    : probeAdaptor->fetch_all_by_sequence('AGTC')
  Description: Fetches all probes with the given sequence
  Returntype : ArrayRef of Bio::EnsEMBL::Funcgen::Probe objects
  Caller     : General
  Status     : Stable

=cut

sub fetch_all_by_sequence {
  my $self = shift;
  my $sequence = shift;

  my $probe_sequence = $self->db->get_ProbeSequenceAdaptor->fetch_by_sequence($sequence);
  return $self->fetch_all_by_ProbeSequence($probe_sequence);
}


=head2 fetch_all_by_ProbeSequence

  Arg [1]    : Bio::EnsEMBL::Funcgen::ProbeSequence
  Example    : probeAdaptor->fetch_all_by_ProbeSequence($probeSequence)
  Description: Fetches all Probes linked to this ProbeSequence object
  Returntype : ArrayRef of Bio::EnsEMBL::Funcgen::Probe objects
  Caller     : General
  Status     : Stable

=cut

sub fetch_all_by_ProbeSequence {
  my ($self, $probe_sequence) = @_;

  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::ProbeSequence', $probe_sequence);
  
  return $self->fetch_all_by_probe_sequence_id($probe_sequence->dbID);
}


=head2 fetch_all_by_probe_sequence_id

  Arg [1]    : String - ProbeSequence dbID
  Example    : probeAdaptor->fetch_all_by_probe_sequence_id
  Description: Fetches all Probes linked to this ProbeSequence dbID
  Returntype : ArrayRef of Bio::EnsEMBL::Funcgen::Probe objects 
  Caller     : General
  Status     : Stable

=cut

sub fetch_all_by_probe_sequence_id {
  my ($self, $probe_sequence_id) = @_;

  $self->bind_param_generic_fetch($probe_sequence_id, SQL_INTEGER);
  return $self->generic_fetch('p.probe_seq_id=?');
}

=head2 fetch_all_by_transcript_stable_id

  Arg [1]    : string - transcript stable id
  Example    : my $probe_list = $probe_adaptor->fetch_all_by_transcript_stable_id('ENST00000489935');
  Description: Fetches all probes that have been mapped to this transcript by the 
               probe2transcript step in the probemapping pipeline.
  Returntype : Arrayref
  Caller     : General

=cut

sub fetch_all_by_transcript_stable_id {
  my $self = shift;
  my $transcript_stable_id = shift;

  my $probe_transcript_mappings = $self->db->get_ProbeTranscriptMappingAdaptor->fetch_all_by_transcript_stable_id($transcript_stable_id);
  
  if (! defined $probe_transcript_mappings) {
    return [];
  }
  
  my @probes_mapped_to_transcript;
  foreach my $current_probe_transcript_mapping (@$probe_transcript_mappings) {
    push @probes_mapped_to_transcript,
      $self->fetch_by_dbID($current_probe_transcript_mapping->probe_id);
  }
  return \@probes_mapped_to_transcript;
}

=head2 fetch_all_by_name

  Arg [1]    : string - probe name
  Example    : my @probes = @{$opa->fetch_all_by_name('Probe1')};
  Description: Convenience method to re-instate the functionality of
               $core_dbentry_adpator->fetch_all_by_External_name('probe_name');
               WARNING: This may not be the probe you are expecting as
               probe names are not unqiue across arrays and vendors.
               These should ideally be validated using the attached array
               information or alternatively use fetch_by_array_probe_probeset_name
               Returns a probe with the given name.
  Returntype : Arrayref
  Exceptions : Throws if name not passed
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_name {
  my ($self, $name) = @_;
  throw('Must provide a probe name argument') if ! defined $name;
  $self->bind_param_generic_fetch($name, SQL_VARCHAR);
  return $self->generic_fetch('p.name=?');
}

=head2 fetch_all_by_ProbeSet

  Arg [1]    : Bio::EnsEMBL::ProbeSet
  Example    : my @probes = @{$opa->fetch_all_by_ProbeSet($pset)};
  Description: Fetch all probes in a particular ProbeSet.
  Returntype : Listref of Bio::EnsEMBL::Probe objects
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_ProbeSet {
  my ($self, $probeset) = @_;
  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::ProbeSet', $probeset);
  return $self->generic_fetch('p.probe_set_id = '.$probeset->dbID);
}

=head2 fetch_all_by_Array

  Arg [1]    : Bio::EnsEMBL::Funcgen::Array
  Example    : my @probes = @{$opa->fetch_all_by_Array($array)};
  Description: Fetch all probes on a particular array.
  Returntype : Listref of Bio::EnsEMBL::Probe objects.
  Exceptions : throws if arg is not valid or stored
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_Array {
  my ($self,$array) = @_;

  $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::Array', $array);

#  if(! (ref($array) && $array->isa('Bio::EnsEMBL::Funcgen::Array') && $array->dbID())) {
#    throw('Need to pass a valid stored Bio::EnsEMBL::Funcgen::Array');
#  }

  return $self->generic_fetch('p.array_id = '.$array->dbID);
}

=head2 fetch_all_by_ArrayChip

  Arg [1]    : Bio::EnsEMBL::Funcgen::ArrayChip
  Example    : my @probes = @{$opa->fetch_all_by_ArrayChip($array_chip)};
  Description: Fetch all probes on a particular ArrayChip.
  Returntype : Listref of Bio::EnsEMBL::Probe objects.
  Exceptions : throw is arg is not valid or stored
  Caller     : General
  Status     : Deprecated

=cut

sub fetch_all_by_ArrayChip {
  my ($self, $array_chip) = @_;
  deprecate('Will be removed in e94. Probes are now only part of 1 array, use $pa->fetch_all_by_Array');

  if(! (ref($array_chip) && $array_chip->isa('Bio::EnsEMBL::Funcgen::ArrayChip') && $array_chip->dbID())){
    throw('Need to pass a valid stored Bio::EnsEMBL::Funcgen::ArrayChip');
  }

  return $self->generic_fetch("p.array_chip_id =".$array_chip->dbID());
}


=head2 fetch_by_ProbeFeature

  Arg [1]    : Bio::EnsEMBL::Funcgen::ProbeFeature
  Example    : my $probe = $opa->fetch_by_ProbeFeature($feature);
  Description: Returns the probe that created a particular feature.
  Returntype : Bio::EnsEMBL::Probe
  Exceptions : Throws if argument is not a Bio::EnsEMBL::Funcgen::ProbeFeature object
  Caller     : General
  Status     : At Risk

=cut

sub fetch_by_ProbeFeature {
  my $self    = shift;
  my $feature = shift;

  if (! ref($feature) ||
      ! $feature->isa('Bio::EnsEMBL::Funcgen::ProbeFeature') ||
      ! $feature->{'probe_id'}
     ) {
    throw('fetch_by_ProbeFeature requires a stored Bio::EnsEMBL::Funcgen::ProbeFeature object');
  }

  return $self->fetch_by_dbID($feature->{'probe_id'});
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
  return (['probe', 'p']);
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
      p.probe_id 
      p.probe_set_id 
      p.name 
      p.length 
      p.array_id 
      p.class 
      p.description 
      p.probe_seq_id
      );
}

=head2 _objs_from_sth

  Arg [1]    : DBI statement handle object
  Example    : None
  Description: PROTECTED implementation of superclass abstract method.
               Creates Probe objects from an executed DBI statement
               handle.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::Probe objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _objs_from_sth {
  my ($self, $sth) = @_;

  my (@result, $array_id, $probe_id, $probe_set_id, $name, $class, $probelength, $desc, $probe_seq_id);
  
  my $array_a = $self->db->get_ArrayAdaptor;
  my $ps_a    = $self->db->get_ProbeSetAdaptor;
  # Caching 
  my $ps_cache    = {};
  my $array_cache = {};
 
  $sth->bind_columns(\$probe_id, \$probe_set_id, \$name, \$probelength, \$array_id, \$class, \$desc, \$probe_seq_id);

  while ( $sth->fetch() ) {

    # Relying on that there is a 1-1 relationship array_chip <-> array
    $array_cache->{$array_id} = 
      $array_a->fetch_by_dbID($array_id) if(!defined $array_cache->{$array_id});
    if(defined $probe_set_id){
      if(!defined $ps_cache->{$probe_set_id}){
        $ps_cache->{$probe_set_id} = $ps_a->fetch_by_dbID($probe_set_id);
      }
    }
    my $array = $array_cache->{$array_id};
    my $probe_set = $ps_cache->{$probe_set_id};
    
    my $probe = Bio::EnsEMBL::Funcgen::Probe->new(
        -dbID          => $probe_id,
        -name          => $name,
        -array         => $array,
        -probe_set     => $probe_set,
        -length        => $probelength,
        -class         => $class,
        -description   => $desc,
        -probe_seq_id  => $probe_seq_id,
        -adaptor       => $self,
        );
    push @result, $probe;
  }
  return \@result;
}

=head2 store

  Arg [1]    : List of Bio::EnsEMBL::Funcgen::Probe objects
  Example    : $opa->store($probe1, $probe2, $probe3);
  Description: Stores given Probe objects in the database. Should only be
               called once per probe because no checks are made for duplicates
                           Sets dbID and adaptor on the objects that it stores.
  Returntype : ARRAYREF
  Exceptions : Throws if arguments are not Probe objects
  Caller     : General
  Status     : At Risk

=cut

sub store {
  my ($self, @probes) = @_;
  
  throw('Must call store with a list of Probe objects') if (scalar @probes == 0);

  my $db = $self->db();
  my $new_sth = $self->prepare
    (
     "INSERT INTO probe( probe_set_id, name, length, array_id, class, description, probe_seq_id)".
     "VALUES (?, ?, ?, ?, ?, ?, ?)"
    );
  my $probe_seq_sth = $self->prepare('INSERT INTO probe_seq (probe_sha1, probe_dna) VALUES (CAST(sha1(?) AS CHAR), ?)');
  my $probe_sequence_adaptor = $self->db->get_ProbeSequenceAdaptor;

 PROBE: foreach my $probe (@probes) {
     
    if ( !ref $probe || ! $probe->isa('Bio::EnsEMBL::Funcgen::Probe') ) {
      throw("Probe must be an Probe object ($probe)");
    }
    
    $self->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::Array', $probe->Array);

    if ( $probe->is_stored($db) ) {
      warning('Probe [' . $probe->dbID() . '] is already stored in the database');
      next PROBE;
    }
# ToDo: Storing a new Probe I don't have an adaptor. I assume a lot of this is coming from the Pipeline?
    my $probe_sequence = '';
    if(! defined $probe->{'probe_sequence'}){
      if(defined $probe->_probe_seq_id){
        say ref($probe_sequence_adaptor);
        $probe_sequence = $probe_sequence_adaptor->fetch_by_dbID($probe->_probe_seq_id);
        say "Done";
      }
      else{
        throw('No probe sequence or ID. This should not be possible');
      }
    }

# ToDo: Replaced by above    
#    my $probe_sequence = $probe->get_ProbeSequence;
    
    if (! defined $probe_sequence->sequence) {
      throw( "Probe sequence not defined in probe:\n" . Dumper($probe));
    }
    
    $probe_sequence_adaptor->store($probe_sequence);
    my $probe_seq_id = $probe_sequence->dbID;
    
# ToDo: Can this happen?
    if (! defined $probe_seq_id) {
      throw( "No probe sequence id for probe::\n" . Dumper($probe));
    }
    $probe->set_ProbeSequence($probe_sequence);

    my $ps_id = (defined $probe->probe_set()) ? $probe->probe_set()->dbID() : undef;
    $new_sth->bind_param(1, $ps_id,              SQL_INTEGER);
    $new_sth->bind_param(2, $probe->name,        SQL_VARCHAR);
    $new_sth->bind_param(3, $probe->length(),    SQL_INTEGER);
    $new_sth->bind_param(4, $probe->Array->dbID, SQL_INTEGER);
    $new_sth->bind_param(5, $probe->class(),     SQL_VARCHAR);
    $new_sth->bind_param(6, $probe->description, SQL_VARCHAR);
    $new_sth->bind_param(7, $probe_seq_id,       SQL_INTEGER);
    $new_sth->execute();
    $probe->dbID($self->last_insert_id);
    $probe->adaptor($self);
  }
  
  return \@probes;
}

1;

