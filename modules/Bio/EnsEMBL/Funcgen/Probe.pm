#
# Ensembl module for Bio::EnsEMBL::Funcgen::Probe
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

Bio::EnsEMBL::Funcgen::Probe - A module to represent a nucleotide probe.

=head1 SYNOPSIS

  use Bio::EnsEMBL::Funcgen::Probe;

  my $probe = Bio::EnsEMBL::Funcgen::Probe->new(
    -PROBE_SET     => $probe_set,
    -NAME          => 'Probe-1',
    -ARRAY         => $array,
    -ARRAY_CHIP_ID => $ac_dbid,
    -CLASS         => "EXPERIMENTAL",
  );

=head1 DESCRIPTION

An Probe object represents an probe on a microarray. The data (currently the
name, probe_set_id, length, pair_index and class) are stored
in the oligo_probe table.

For Affy arrays, a probe can be part of more than one array, but only part of
one probeset. On each Affy array the probe has a slightly different name. For
example, two different complete names for the same probe might be
DrosGenome1:AFFX-LysX-5_at:535:35; and Drosophila_2:AFFX-LysX-5_at:460:51;. In
the database, these two probes will have the same oligo_probe_id. Thus the same
Affy probe can have a number of different names and complete names depending on
which array it is on.

=cut

package Bio::EnsEMBL::Funcgen::Probe;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument qw( rearrange ) ;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Utils::Exception qw( deprecate );
use Bio::EnsEMBL::Funcgen::ProbeSequence;

use feature qw(say);
use Data::Dumper;

use base qw( Bio::EnsEMBL::Funcgen::Storable );


=head2 new

  Arg [-NAME]          : string - probe name
        Used when the probe is on one array.
  Arg [-ARRAY]          : Bio::EnsEMBL::Funcgen::Array
        Used when the probe is on one array.
  Arg [-NAMES]          : Listref of ints - array_chip db IDs
        Used when the probe is on multiple arrays.
  Arg [-PROBE_SET]      : Bio::EnsEMBL::ProbeSet
        Each probe is part of one(and only one) probeset, if not probe set
        then probeset = probe i.e. probe_set size = 1
  Arg [-LENGTH]         : int - probe length
        Will obviously be the same for all probes if same probe
		is on multiple arrays.
  Arg [-CLASS]          : string - probe class e.g. CONTROL, EXPERIMENTAL
        Will be the same for all probes if same probe is on
		multiple arrays.
  Arg [-DESCRIPTION]    : (optional) string - description


  Example    : my $probe = Bio::EnsEMBL::Probe->new(
      -NAME          => 'Probe-1',
      -PROBE_SET     => $probe_set,
      -ARRAY         => $array,
      -ARRAY_CHIP_ID => $array_chip_id,
      -LENGTH        => 25,
      -CLASS         => 'EXPERIMENTAL',
      -DESCRIPTION   => 'Some useful description',
    );
  Description: Creates a new Bio::EnsEMBL::Probe object.
  Returntype : Bio::EnsEMBL::Probe
  Exceptions : Throws if not supplied with probe name and array
  Caller     : General
  Status     : Medium Risk

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;

  my $self = $class->SUPER::new(@_);
  

  my (
      $name,
      $array,
      $probe_set,       
      $aclass,
      $length,         
      $desc,
      $sequence,
      $probe_seq_id
     ) = rearrange([
      'NAME',
      'ARRAY',
      'PROBE_SET',      
      'CLASS',
      'LENGTH',         
      'DESCRIPTION',
      'SEQUENCE',
      'PROBE_SEQ_ID'
      ], @_);

  if(!defined $array) {
    throw("Bio::EnsEMBL::Funcgen::Array required");
  }
  if(defined($array) && (!ref($array) || !$array->isa('Bio::EnsEMBL::Funcgen::Array'))) {
    throw("Not a Bio::EnsEMBL::Funcgen::Array object");
  }

  $self->{array} = $array;

  $self->name($name);

  $self->probe_set($probe_set)  if defined $probe_set;
  $self->class($aclass)         if defined $aclass;
  $self->length($length)        if defined $length;
  $self->description($desc)     if defined $desc;

  if(defined $probe_seq_id and defined $sequence){
    throw("Either define -SEQUENCE or PROBE_SEQ_ID, not both");
  }

  if(!defined $probe_seq_id and !defined $sequence){
    throw("Define either -SEQUENCE or PROBE_SEQ_ID");
  }

  $self->_probe_seq_id($probe_seq_id) if defined $probe_seq_id;

  if (defined $sequence) {
    my $probe_sequence = Bio::EnsEMBL::Funcgen::ProbeSequence->new(
      -sequence => $sequence
    );
    $self->set_ProbeSequence($probe_sequence);
  }

  return $self;
}

=head2 sequence

  Arg [1]    : Optional - String sequence
  Example    : $probe->sequence
  Description: Getter/Setter for probe sequence
  Returntype : String - sequence
  Caller     : General
  Status     : Stable

=cut

sub sequence {
    my $self = shift;
    $self->{sequence} = shift if @_;

    if (! defined $self->{sequence}) {
      $self->{sequence} = $self->get_ProbeSequence->sequence;
    }

    return $self->{sequence};
}

=head2 

  Arg [1]    : Optional - String _probe_seq_id
  Example    : $self->_probe_seq_id
  Description: Getter/Setter for probe_seq_id
  Returntype : String - probe_seq_id
  Caller     : Private
  Status     : Stable

=cut

sub _probe_seq_id {
    my $self = shift;
    $self->{probe_seq_id} = shift if @_;
    return $self->{probe_seq_id};
}

=head2 name

  Arg [1]    : Optional - String - name of probe
  Example    : $probe->name
  Description: Getter/Setter for name
  Returntype : String - Name of the probe
  Caller     : General
  Status     : Stable

=cut

sub name {
    my $self = shift;
    $self->{name} = shift if @_;
    return $self->{name};
}


=head2 

  Arg [1]    : None
  Example    : $probe->get_ProbeSequence
  Description: Getter for Probe->ProbeSequence
  Returntype : String - sequence
  Caller     : General
  Status     : Stable

=cut

sub get_ProbeSequence {
    my $self = shift;

    if (
         (! defined $self->{'probe_sequence'})
      && (defined $self->_probe_seq_id)
    ) {
        $self->{'probe_sequence'} = $self->_fetch_ProbeSequence;
    }
    return $self->{'probe_sequence'};
}

=head2 _fetch_ProbeSequence

  Arg [1]    : None
  Example    : $self->_fetch_ProbeSequence
  Description: Fetches ProbeSequence linked to this Probe
  Returntype : String - probe sequence
  Caller     : Private
  Status     : Stable

=cut

sub _fetch_ProbeSequence {
    my $self = shift;
    return $self->adaptor()->db()->get_ProbeSequenceAdaptor()->fetch_by_dbID($self->_probe_seq_id);
}

=head2 set_ProbeSequence

  Arg [1]    : String sequence
  Example    : $probe->set_ProbeSequence('AGTC')
  Description: Sets sequence for this Probe
  Returntype : None
  Caller     : General
  Status     : Stable

=cut

sub set_ProbeSequence {
    my ($self, $probe_sequence) = @_;
    throw 'No sequence provided' if(! defined $probe_sequence);
    $self->{'probe_sequence'} = $probe_sequence;
}



=head2 new_fast

  Args       : Hashref with all internal attributes set
  Example    : none
  Description: Quick and dirty version of new. Only works if the code is very
               disciplined. Cannot add array chip probe names unless we recreate
               the data structure in the caller.
  Returntype : Bio::EnsEMBL::Funcgen::Probe
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub new_fast {
  bless ($_[1], $_[0]);
}


=head2 get_all_ProbeFeatures

  Args       : None
  Example    : my $features = $probe->get_all_ProbeFeatures();
  Description: Get all features produced by this probe. The probe needs to be
               database persistent.
  Returntype : Listref of Bio::EnsEMBL:Funcgen::ProbeFeature objects
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub get_all_ProbeFeatures {
  my $self = shift;
  if ( $self->adaptor() && $self->dbID() ) {
    return $self->adaptor()->db()->get_ProbeFeatureAdaptor()->fetch_all_by_Probe($self);
  } else {
    warning('Need database connection to retrieve Features');
    return [];
  }
}

=head2 Array

  Args       : None
  Example    : my $array = $probe->Array();
  Description: Returns the array that this probe is part of. Only works if the
               probe was retrieved from the database or created using
               add_Array_probename (rather than add_arrayname_probename).
  Returntype : Bio::EnsEMBL::Funcgen::Array 
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub Array {
  my $self = shift;

 return $self->{array};
}

####### Boulevard of Broken Dreams (Deprecated methods)

=head2 _array_chip_id

  Arg [1]    : Optional
  Example    : $probe->
  Description: 
  Returntype : 
  Caller     : General
  Status     : Deprecate

=cut

sub _array_chip_id {
    my $self = shift;
    deprecate('Probes are linked to one Array only. array_chip is not needed anymore');
    $self->{'_array_chip_id'} = shift if @_;
    return $self->{'_array_chip_id'};
}

=head2 array_chip

  Arg [1]    : Optional - ArrayChip dbID
  Example    : $probe->
  Description: Getter/Setter for ArrayChip dbID for this Probe
  Returntype : String - ArrayChip dbID
  Caller     : General
  Status     : Stable

=cut

sub array_chip {
    my ($self, $array_chip) = @_;

    if (defined $array_chip) {
      $self->{_array_chip} = $array_chip;
    }
    if (! defined $self->{_array_chip}) {
      $self->{_array_chip} = $self->adaptor->db->get_ArrayChipAdaptor->fetch_by_dbID($self->_array_chip_id);
    }
    if (! defined $self->{_array_chip}) {
      die("Probe was not linked to an array chip!");
    }

    return $self->{_array_chip};
}

=head2 add_array_chip_probename

  Arg [1]    : string - probe name
  Arg [2]    : Bio::EnsEMBL::Funcgen::Array
  Example    : $probe->add_array_chip_probename($probename, $array);
  Description: Adds a probe name / array pair to a probe, allowing incremental
               generation of a probe.
  Returntype : None
  Exceptions : None
  Caller     : General,
               Probe->new(),
               ProbeAdaptor->_obj_from_sth(),
               AffyProbeAdaptor->_obj_from_sth()
  Status     : Deprecated

=cut

sub add_array_chip_probename {
  deprecate('Will be removed in e94. Probes are linked to one Array only. array_chip is not needed anymore');
  my $self = shift;
  my ($probename, $array) = @_;
  $self->{arrays}     ||= {};
  $self->{probenames} ||= {};

  if(! (ref($array) && $array->isa('Bio::EnsEMBL::Funcgen::Array'))){
    throw('You must pass a valid Bio::EnsEMBL::Funcgen::Array. ')
  }
  
  $self->{arrays}->{$array->name}           = $array;
#   $self->{arrays}->{$ac_dbid}           = $array;
  $self->{probenames}->{$array->name} ||= [];
  push @{$self->{probenames}->{$array->name}}, $probename;
#say Dumper($self);
  return;
#  $self->{arrays}     ||= {};
#  $self->{probenames} ||= {};
#
#  if(! (ref($array) && $array->isa('Bio::EnsEMBL::Funcgen::Array'))){
#    throw('You must pass a valid Bio::EnsEMBL::Funcgen::Array. ')
#  }
#
#  $self->{arrays}->{$array->name}           = $array;
#  $self->{probenames}->{$array->name} ||= [];
#  push @{$self->{probenames}->{$array->name}}, $probename;
#
#  return;
}


=head2 get_all_Arrays

  Args       : None
  Example    : my $arrays = $probe->get_all_Arrays();
  Description: Returns all arrays that this probe is part of. Only works if the
               probe was retrieved from the database or created using
               add_Array_probename (rather than add_arrayname_probename).
  Returntype : Listref of Bio::EnsEMBL::Funcgen::Array objects
  Exceptions : None
  Caller     : General
  Status     : Depracted

=cut

sub get_all_Arrays {
  my $self = shift;
  deprecate('Will be removed in e94. Probes are now only part of 1 array, use $probe->array');
  return [$self->{array}];
  
}


=head2 get_names_Arrays

  Args       : None
  Example    : my %name_array_pairs = %{$probe->get_names_Arrays};
  Description: Returns Array name hash
  Returntype : hashref of probe name Bio::EnsEMBL::Funcgen::Array pairs
  Exceptions : None
  Caller     : General
  Status     : Deprecated

=cut

sub get_names_Arrays {
  my $self = shift;
  deprecate('Will be removed in e94. Probes are now only part of 1 array, use $probe->Array->name instead ');
  return {$self->{array}->name,$self->{array}};
}

=head2 get_all_probenames

  Arg [1]    : Optional - list of array names, defaults to all available
  Example    : my @probenames = @{$probe->get_all_probenames()};
  Description: Retrieves all names for this probe. Only makes sense for probes
               which share identical sequence for a given probeset and array.
               This can either be:
               1 A non-probeset array where probes with different names but identical
                 sequence have beem merged, this is only possible if the probes in question
                 share the same Array but are on seaprate ArrayChips.
               2 If they are part of a probeset (i.e. Affy probes), in which case
                 get_all_complete_names() would be more appropriate.
  Returntype : Arrayref of strings
  Exceptions : None
  Caller     : General
  Status     : Deprecated

=cut

sub get_all_probenames {
  my $self        = shift;
  deprecate('Will be removed in e94. Probes are unique within an array');
  return([$self->name]);
  my @array_names = @_;

  if (! @array_names) {
    @array_names = keys %{$self->{'probenames'}};
  }

  my @probe_names;
  foreach my $current_array_name (@array_names) {
    push @probe_names, @{$self->{'probenames'}->{$current_array_name}};
  }
  return \@probe_names;
}

=head2 get_probename

  Arg [1]    : string - array name
  Example    : my $probename = $probe->get_probename('Array-1');
  Description: For a given array, retrieve the name for this probe.
  Returntype : string
  Exceptions : Throws if the array name is required but not specified
               Warns if probe has more than one name for the given array.
  Caller     : General
  Status     : Deprecate

=cut

sub get_probename {
  my ($self, $arrayname) = @_;
  deprecate('Will be removed in e94. Probes are unique within an array, use $probe->name instead');
  return $self->{name};
}



=head2 get_all_complete_names

  Args       : None
  Example    : my @compnames = @{$probe->get_all_complete_names()};
  Description: Retrieves all complete names for this probe. The complete name
               is a concatenation of the array name, the probeset name and the
               probe name.
  Returntype : Arrayref of strings
  Exceptions : None
  Caller     : Used by web for the names like here: http://www.ensembl.org/Homo_sapiens/Transcript/Oligos?db=core;g=ENSG00000139618;r=13:32315474-32400266;t=ENST00000470094
  Status     : Deprecated

=cut

sub get_all_complete_names {
  my $self = shift;

  deprecate('Will be removed in e94. Probes are unique now in an Array and ProbeSet. Use $probe->get_full_name instead');
  my $probeset;

  if (defined $self->probeset && $self->probeset->name) {
    $probeset = ':' . $self->probeset->name . ':';
  } else {
    $probeset = ':';
  }

  my $complete_name = $self->Array->name . $probeset . $self->name;
  return [$complete_name];
}

=head2 get_full_name

  Args       : None
  Example    : my $full_name = $probe->get_full_name;
  Description: Retrieves the full name for this probe. The complete name
               is a concatenation of the array name, the probeset name and the
               probe name.
  Returntype : string
  Exceptions : None
  Caller     : Used by web for the names like here: http://www.ensembl.org/Homo_sapiens/Transcript/Oligos?db=core;g=ENSG00000139618;r=13:32315474-32400266;t=ENST00000470094
  Status     : Stable

=cut

sub get_full_name {
  my $self = shift;

  my $probeset;

  if (defined $self->probeset && $self->probeset->name) {
    $probeset = ':' . $self->probeset->name . ':';
  } else {
    $probeset = ':';
  }

  my $full_name = $self->Array->name . $probeset . $self->name;
  return $full_name;
}

=head2 get_complete_name

   Arg [1]    : string - array name
   Example    : my $compname = $probe->get_complete_name('Array-1');
   Description: For a given array, retrieve the complete name for this probe.
   Returntype : string
   Exceptions : Throws if the array name not specified or not known for this probe
   Caller     : General
   Status     : Deprecated

=cut

sub get_complete_name {
   my $self = shift;
   my $arrayname = shift;

    deprecate(
        "get_complete_name has been deprecated and will be removed in Ensembl
        release 93."
    );

   throw('Must provide and array name argument to retreive the complete name') if ! defined $arrayname;

   my $probename = $self->get_probename($arrayname);

   if (!defined $probename) {
     throw('Unknown array name');
   }

   my $probeset = $self->probeset()->name();
   $probeset .= ':' if $probeset;

   return "$arrayname:$probeset$probename";
}

=head2 probe_set

  Arg [1]    : (optional) Bio::EnsEMBL::Funcgen::ProbeSet
  Example    : my $probe_set = $probe->probeset();
  Description: Getter and setter of probe_set attribute for Probe objects.
  Returntype : Bio::EnsEMBL::Funcgen::ProbeSet
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub probe_set {
    my $self = shift;

    $self->{'probe_set'} = shift if @_;
    return $self->{'probe_set'};
}

sub probeset {
    my $self = shift;
    return $self->probe_set(@_);
}

=head2 class

  Arg [1]    : (optional) string - class
  Example    : my $class = $probe->class();
  Description: Getter and setter of class attribute for Probe
               objects e.g. CONTROL, EXPERIMENTAL
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub class {
    my $self = shift;
    $self->{'class'} = shift if @_;
    return $self->{'class'};
}

=head2 length

  Arg [1]    : (optional) int - probe length
  Example    : my $probelength = $probe->length();
  Description: Getter and setter of length attribute for Probe
               objects.
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub length {
  my $self = shift;
  $self->{'length'} = shift if @_;

  if (! defined $self->{'length'}) {
    $self->{'length'} = length($self->sequence);
  }
  return $self->{'length'};
}

=head2 description

  Arg [1]    : (optional) string - description
  Example    : my $pdesc = $probe->description();
  Description: Getter and setter of description attribute for Probe
               objects.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub description {
  my $self = shift;
  $self->{'description'} = shift if @_;
  return $self->{'description'};
}


# =head2 feature_count
#
#   Arg[0]     : recount flag
#   Example    : my $num_features = $probe->feature_count();
#   Description: Counts the number of ProbeFeatures associated with this Probe
#   Returntype : int
#   Exceptions : None
#   Caller     : General
#   Status     : Medium Risk
#
# =cut
#
#
# sub feature_count{
#   my ($self, $recount) = @_;
#
#   if($recount ||
#     (! $self->{feature_count})){
#     $self->{feature_count} = $self->adaptor->db->get_ProbeFeatureAdaptor->count_probe_features_by_probe_id($self->dbID);
#   }
#
#   return $self->{feature_count};
# }

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

sub get_all_Transcript_DBEntries {
  my $self = shift;
  deprecate(
    "get_all_Transcript_DBEntries has been deprecated and will be removed in Ensembl release 92."
        . " Please use fetch_all_ProbeTranscriptMappings instead."
  );
  return $self->fetch_all_ProbeTranscriptMappings;
}

=head2 fetch_all_ProbeTranscriptMappings

  Arg[0]     : none
  Example    : $probe->fetch_all_mapped_Transcripts;
  Description: Returns all mappings of this probe to transcripts.
  Returntype : Listref of Bio::EnsEMBL::Funcgen::ProbeTranscriptMapping objects
  Exceptions : none
  Caller     : general
  Status     : at risk

=cut

sub fetch_all_ProbeTranscriptMappings {
  my $self = shift;
  return $self->adaptor->db->get_ProbeTranscriptMappingAdaptor->fetch_all_by_probe_id($self->dbID);
}

1;
