#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::ProbeAdaptor
#

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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

use base qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor);

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
  Example    : 
                my $probe_set_array_name          = 'HG-U133A';
                my $probe_name_on_probe_set_array = '552:605;';
                my $probe_set_name                = '214727_at';

                my $probe
                  = $probe_adaptor->fetch_by_array_probe_probeset_name(
                  
                    $probe_set_array_name, 
                    $probe_name_on_probe_set_array, 
                    $probe_set_name
                  );

  Description: Returns a probe given a combination of array name, probeset and
               probe name. This will uniquely define an Affy probe. Only one 
               probe is ever returned.
  Returntype : Bio::EnsEMBL::Funcgen::Probe
  Exceptions : throws if array or probe name not defined
  Caller     : General
  Status     : At Risk - rename to fetch_by_probe_array_probeset_name?

=cut
sub fetch_by_array_probe_probeset_name {

  my $self          = shift;
  my $array_name    = shift;
  my $probe_name    = shift;
  my $probeset_name = shift;

  if ($probeset_name) {
    return 
      $self->fetch_by_array_probe_set_probe_name(
        $array_name, 
        $probeset_name, 
        $probe_name
      );
  }
  return 
    $self->fetch_by_array_probe_name(
      $array_name, 
      $probe_name
    );
}

=head2 fetch_by_array_probe_name

  Arg [1]    : string - name of array
  Arg [2]    : string - name of probe

  Example    : 
                my $array_name = 'OneArray';
                my $probe_name = 'PH_hs_0000014';

                my $probe = $probe_adaptor->fetch_by_array_probe_name(
                  $array_name,
                  $probe_name
                );

  Description: Returns a probe given a combination of array and probe 
               name. Use this to fetch probes from arrays that are not
               organized as probe sets.
               
  Returntype : Bio::EnsEMBL::Funcgen::Probe
  Caller     : General

=cut
sub fetch_by_array_probe_name {

  my $self       = shift;
  my $array_name = shift;
  my $probe_name = shift;
  
  my $sql = "
    select 
      distinct(probe.probe_id) 
    from 
      probe 
      join array_chip using (array_chip_id)
      join array      using (array_id)
    where
      array.name         = ?
      and probe.name     = ?
  ";
  
  my $dbc = $self->db->dbc;
  my $sth = $dbc->prepare($sql);
  
  $sth->bind_param(1, $array_name,    SQL_VARCHAR);
  $sth->bind_param(2, $probe_name,    SQL_VARCHAR);

  $sth->execute;
  
  my @probe_ids = $sth->fetchrow_array;
  my $number_of_probes_found = @probe_ids;
  my ($probe_id) = @probe_ids;
  
  if ($number_of_probes_found == 0) {
    return;
  }
  if ($number_of_probes_found > 1) {
    warn "Found $number_of_probes_found probes with name $probe_name on array ${array_name}. Returning the first one found.";
  }
  my $probe = $self->fetch_by_dbID($probe_id);
  return $probe;
}

=head2 fetch_by_array_probe_set_probe_name

  Arg [1]    : string - name of array
  Arg [2]    : string - name of probe set
  Arg [3]    : string - name of probe

  Example    : 
                my $probe_set_array_name          = 'HuGene-2_0-st-v1';
                my $probe_name_on_probe_set_array = '2398055-16657436';
                my $probe_set_name                = '16657436';

                my $probe
                  = $probe_adaptor->fetch_by_array_probe_set_probe_name(

                      $probe_set_array_name, 
                      $probe_set_name,
                      $probe_name_on_probe_set_array, 
                  );

  Description: Returns a probe given a combination of array and probe 
               name. Use this to fetch probes from arrays that are organized 
               as probe sets.

  Returntype : Bio::EnsEMBL::Funcgen::Probe
  Caller     : General

=cut
sub fetch_by_array_probe_set_probe_name {

  my $self           = shift;
  my $array_name     = shift;
  my $probe_set_name = shift;
  my $probe_name     = shift;

  my $sql = "
    select 
      distinct(probe.probe_id) 
    from 
      probe 
      join array_chip using (array_chip_id)
      join array      using (array_id)
      join probe_set  using (probe_set_id)
    where
      array.name         = ?
      and probe.name     = ?
      and probe_set.name = ?
  ";
  
  my $dbc = $self->db->dbc;
  my $sth = $dbc->prepare($sql);
  
  $sth->bind_param(1, $array_name,     SQL_VARCHAR);
  $sth->bind_param(2, $probe_name,     SQL_VARCHAR);
  $sth->bind_param(3, $probe_set_name, SQL_VARCHAR);

  $sth->execute;
  
  my @probe_ids = $sth->fetchrow_array;
  my $number_of_probes_found = @probe_ids;
  my ($probe_id) = @probe_ids;
  
  if ($number_of_probes_found == 0) {
    return;
  }
  if ($number_of_probes_found > 1) {
    warn "Found $number_of_probes_found probes with name $probe_name on array ${array_name}. Returning the first one found.";
  }
  my $probe = $self->fetch_by_dbID($probe_id);
  return $probe;
}

sub fetch_all_by_sequence {
  my $self = shift;
  my $sequence = shift;

  my $probe_sequence = $self->db->get_ProbeSequenceAdaptor->fetch_by_sequence($sequence);
  return $self->fetch_all_by_ProbeSequence($probe_sequence);
}

sub fetch_all_by_ProbeSequence {
  my $self = shift;
  my $probe_sequence = shift;
  
  if (! defined $probe_sequence->dbID) {
    die;
  }
  return $self->fetch_all_by_probe_sequence_id($probe_sequence->dbID);
}

sub fetch_all_by_probe_sequence_id {
  my $self = shift;
  my $probe_sequence_id = shift;

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
  my $self  = shift;
  my $array = shift;

  if(! (ref($array) && $array->isa('Bio::EnsEMBL::Funcgen::Array') && $array->dbID())) {
    throw('Need to pass a valid stored Bio::EnsEMBL::Funcgen::Array');
  }

  return $self->generic_fetch('p.array_chip_id IN ('.join(',', @{$array->get_array_chip_ids()}).')');
}

=head2 fetch_all_by_ArrayChip

  Arg [1]    : Bio::EnsEMBL::Funcgen::ArrayChip
  Example    : my @probes = @{$opa->fetch_all_by_ArrayChip($array_chip)};
  Description: Fetch all probes on a particular ArrayChip.
  Returntype : Listref of Bio::EnsEMBL::Probe objects.
  Exceptions : throw is arg is not valid or stored
  Caller     : General
  Status     : At Risk

=cut

sub fetch_all_by_ArrayChip {
  my $self  = shift;
  my $array_chip = shift;

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
  return qw( p.probe_id p.probe_set_id p.name p.length p.array_chip_id p.class p.description p.probe_seq_id);
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

  my (@result, $current_dbid, $arraychip_id, $probe_id, $probe_set_id, $name, $class, $probelength, $desc, $probe_seq_id);
  my $array;
  
  my $array_cache     = $self->{_array_cache};
  my $probe_set_cache = $self->{_probe_set_cache};
  
  $sth->bind_columns(\$probe_id, \$probe_set_id, \$name, \$probelength, \$arraychip_id, \$class, \$desc, \$probe_seq_id);

  my $probe;
  while ( $sth->fetch() ) {

    if (! exists $array_cache->{$arraychip_id}) {
      $array_cache->{$arraychip_id} = $self->db->get_ArrayAdaptor()->fetch_by_array_chip_dbID($arraychip_id);
    }
    $array = $array_cache->{$arraychip_id};
    
    my $probe_set;

    if($probe_set_id) {
    
      if (! exists $probe_set_cache->{$probe_set_id}) {
        $probe_set_cache->{$probe_set_id} = $self->db->get_ProbeSetAdaptor()->fetch_by_dbID($probe_set_id);
      }
      $probe_set = $probe_set_cache->{$probe_set_id};
    }

    if (!$current_dbid || $current_dbid != $probe_id) {

      $probe = Bio::EnsEMBL::Funcgen::Probe->new(
        -dbID          => $probe_id,
        -name          => $name,
        -array_chip_id => $arraychip_id,
        -array         => $array,
        -probe_set     => $probe_set,
        -length        => $probelength,
        -class         => $class,
        -description   => $desc,
        -probe_seq_id  => $probe_seq_id,
        -adaptor       => $self,
      );
      push @result, $probe;
      $current_dbid = $probe_id;
    } else {
      $probe->add_array_chip_probename( $name, $array);
    }
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
  my $db = $self->db();
  throw('Must call store with a list of Probe objects') if (scalar @probes == 0);

  my $new_sth = $self->prepare
    (
     "INSERT INTO probe( probe_set_id, name, length, array_chip_id, class, description, probe_seq_id)".
     "VALUES (?, ?, ?, ?, ?, ?, ?)"
    );
  
  my $existing_sth = $self->prepare
    (
     "INSERT INTO probe( probe_id, probe_set_id, name, length, array_chip_id, class, description, probe_seq_id )".
     "VALUES (?, ?, ?, ?, ?, ?, ?, ?)"
    );

  my $probe_seq_sth = $self->prepare('insert into probe_seq (probe_sha1, probe_dna) values (cast(sha1(?) as char), ?)');
  
#   my $array_chip_adaptor = $self->db->get_ArrayChipAdaptor;
#   print Dumper($array_chip_adaptor);
  my $probe_sequence_adaptor = $self->db->get_ProbeSequenceAdaptor;

 PROBE: foreach my $probe (@probes) {
     
    if ( !ref $probe || ! $probe->isa('Bio::EnsEMBL::Funcgen::Probe') ) {
      throw("Probe must be an Probe object ($probe)");
    }
    
    if ( $probe->is_stored($db) ) {
      warning('Probe [' . $probe->dbID() . '] is already stored in the database');
      next PROBE;
    }
    
    # ------------------------------------------------------------------------------------------
    
    my $probe_sequence = $probe->get_ProbeSequence;
    
    if (! defined $probe_sequence->sequence) {
      use Carp;
      confess(
        "Probe sequence not defined in probe:\n"
        . Dumper($probe)
      );
    }
    
    $probe_sequence_adaptor->store($probe_sequence);
    my $probe_seq_id = $probe_sequence->dbID;
    
    if (! defined $probe_seq_id) {
      confess(
        "No probe sequence id for probe::\n"
        . Dumper($probe)
      );
    }
    
    # ------------------------------------------------------------------------------------------

    # Get all the arrays this probe is on and check they're all in the database
    my %array_hashes;
    
    foreach my $ac_id (keys %{$probe->{'arrays'}}) {
      
      if (defined ${$probe->{'arrays'}}{$ac_id}->dbID()) {
        #Will this ever work as generally we're creating from scratch 
        #and direct access to keys above by passes DB fetch
        $array_hashes{$ac_id} = $probe->{'arrays'}{$ac_id};
      }
    }

    throw('Probes need attached arrays to be stored in the database') if ( ! %array_hashes );

    # Insert separate entry (with same probe_id) for each array/array_chip the probe is on

    foreach my $ac_id (keys %array_hashes) {

      my $ps_id = (defined $probe->probe_set()) ? $probe->probe_set()->dbID() : undef;

      foreach my $name (@{$probe->get_all_probenames($array_hashes{$ac_id}->name)}) {

        if ($probe->dbID) {    # Already stored
          $existing_sth->bind_param(1, $probe->dbID,        SQL_INTEGER);
          $existing_sth->bind_param(2, $ps_id,              SQL_INTEGER);
          $existing_sth->bind_param(3, $name,               SQL_VARCHAR);
          $existing_sth->bind_param(4, $probe->length(),    SQL_INTEGER);
          $existing_sth->bind_param(5, $probe->array_chip->dbID,              SQL_INTEGER);
          $existing_sth->bind_param(6, $probe->class(),     SQL_VARCHAR);
          $existing_sth->bind_param(7, $probe->description, SQL_VARCHAR);
          $existing_sth->bind_param(8, $probe_seq_id,       SQL_INTEGER);
          $existing_sth->execute();
        } else {
          # New probe
          $new_sth->bind_param(1, $ps_id,              SQL_INTEGER);
          $new_sth->bind_param(2, $name,               SQL_VARCHAR);
          $new_sth->bind_param(3, $probe->length(),    SQL_INTEGER);
          $new_sth->bind_param(4, $probe->array_chip->dbID,              SQL_INTEGER);
          $new_sth->bind_param(5, $probe->class(),     SQL_VARCHAR);
          $new_sth->bind_param(6, $probe->description, SQL_VARCHAR);
          $new_sth->bind_param(7, $probe_seq_id,       SQL_INTEGER);
          $new_sth->execute();
          $probe->dbID($self->last_insert_id);
          $probe->adaptor($self);
        }
      }
    }
  }
  
  return \@probes;
}

1;

