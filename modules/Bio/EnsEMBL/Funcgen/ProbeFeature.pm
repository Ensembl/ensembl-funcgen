#
# Ensembl module for Bio::EnsEMBL::Funcgen::ProbeFeature
#

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::ProbeFeature - A module to represent an nucleotide probe
genomic mapping.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::ProbeFeature;

my $feature = Bio::EnsEMBL::Funcgen::ProbeFeature->new
 (
	-PROBE         => $probe,
	-MISMATCHCOUNT => 0,
	-SLICE         => $chr_1_slice,
	-START         => 1000000,
	-END           => 1000024,
	-STRAND        => -1,
  -ANALYSIS      => $analysis,
  -CIGAR_STRING  => '1U2M426D2M1m21M',
 );

=head1 DESCRIPTION

An ProbeFeature object represents the genomic placement of an Probe
object. The data are stored in the probe_feature table.

=cut

package Bio::EnsEMBL::Funcgen::ProbeFeature;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument          qw( rearrange );
use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( median );

use base qw( Bio::EnsEMBL::Feature Bio::EnsEMBL::Funcgen::Storable );


=head2 new

  Arg [-PROBE]        : Bio::EnsEMBL::Funcgen::Probe - probe
        A ProbeFeature must have a probe. This probe must already be stored if
		you plan to store the feature.
  Arg [-MISMATCHCOUNT]: int
        Number of mismatches over the length of the probe. 
  Arg [-SLICE]        : Bio::EnsEMBL::Slice
        The slice on which this feature is.
  Arg [-START]        : int
        The start coordinate of this feature relative to the start of the slice
		it is sitting on. Coordinates start at 1 and are inclusive.
  Arg [-END]          : int
        The end coordinate of this feature relative to the start of the slice
		it is sitting on. Coordinates start at 1 and are inclusive.
  Arg [-STRAND]       : int
        The orientation of this feature. Valid values are 1, -1 and 0.
  Arg [-dbID]         : (optional) int
        Internal database ID.
  Arg [-ADAPTOR]      : (optional) Bio::EnsEMBL::DBSQL::BaseAdaptor
        Database adaptor.
  Example    : my $feature = Bio::EnsEMBL::Funcgen::ProbeFeature->new(
				   -PROBE         => $probe,
				   -MISMATCHCOUNT => 0,
				   -SLICE         => $chr_1_slice,
				   -START         => 1_000_000,
				   -END           => 1_000_024,
				   -STRAND        => -1,
				   -ANALYSIS      => $analysis,
                   -CIGARLINE     => '15M2m3d4M', 
                   #Can represent transcript alignment as gapped genomic alignments
                   #D(eletions) representing introns
                   #Lowercase m's showing sequence mismatches
			   ); 
  Description: Constructor for ProbeFeature objects.
  Returntype : Bio::EnsEMBL::Funcgen::ProbeFeature
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub new {
  my $caller = shift;
	
  my $class = ref($caller) || $caller;
	
  my $self = $class->SUPER::new(@_);
	
  my ($probe, $mismatchcount, $pid, $cig_line, $hit_id, $source)
    = rearrange(['PROBE', 'MISMATCHCOUNT', 'PROBE_ID', 'CIGAR_STRING', 'HIT_ID', 'SOURCE'], @_);
    
  $self->{'probe_id'} = $pid if $pid;
  $self->probe($probe) if $probe;
  $self->mismatchcount($mismatchcount)  if defined $mismatchcount;#do not remove until probe mapping pipeline fixed
  $self->cigar_string($cig_line)        if defined $cig_line;
  $self->hit_id($hit_id)                if defined $hit_id;
  $self->source($source)                if defined $source;

  return $self;
}

=head2 new_fast

  Args       : Hashref with all internal attributes set
  Example    : none
  Description: Quick and dirty version of new. Only works if the code is very
               disciplined. 
  Returntype : Bio::EnsEMBL::Funcgen::ProbeFeature
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub new_fast {
  bless ($_[1], $_[0]);
}

=head2 probeset

  Arg [1]    : (optional) string - probeset
  Example    : my $probeset = $feature->probeset();
  Description: Getter and setter for the probeset for this feature. Shortcut
               for $feature->probe->probeset(), which should be used instead.
			   Probeset is not persisted if set with this method.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Medium Risk
             : Use $feature->probe->probeset() because this may be removed

=cut

sub probeset {
  my $self = shift;

  $self->{'probeset'} = shift if @_;

  if (! $self->{'probeset'}) {
    $self->{'probeset'} = $self->probe->probe_set;
  }
  return $self->{'probeset'};
}

sub probe_set_id {
  my $self = shift;
  
  if (! defined $self->{'_probe_set_id'}) {
    my $probe_set = $self->probeset;
    if (defined $probe_set) {
      $self->{'_probe_set_id'} = $probe_set->dbID;
    }
  }
  return $self->{'_probe_set_id'};
}

# probe_set is more logical.
sub probe_set {
  my $self = shift;
  return $self->probeset;
}

=head2 mismatchcount

  Arg [1]    : int - number of mismatches
  Example    : my $mismatches = $feature->mismatchcount();
  Description: Getter and setter for number of mismatches for this feature.
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : High Risk

=cut

sub mismatchcount {
    my $self = shift;
    $self->{'mismatchcount'} = shift if @_;
    return $self->{'mismatchcount'};
}

=head2 hit_id

  Arg [1]    : String - identifier of the sequence on which the match from which the probe feature was created was found
  Example    : my $hit_id = $feature->hit_id;
  Description: Getter and setter for the hit id of this feature.
  Returntype : int
  Exceptions : None
  Caller     : General

=cut
sub hit_id {
    my $self = shift;
    $self->{'hit_id'} = shift if @_;
    return $self->{'hit_id'};
}

=head2 source

  Arg [1]    : Enum[String] - Either 'genomic' or 'transcript'
  Example    : my $source = $feature->source;
  Description: Returns 'genomic', if the match for this probe feature was on genomic sequence, 'transcript', if it was found on a transcript.
  Returntype : Enum[String] - Either 'genomic' or 'transcript'
  Exceptions : None
  Caller     : General

=cut
sub source {
    my $self = shift;
    $self->{'source'} = shift if @_;
    return $self->{'source'};
}

=head2 cigar_string

  Arg [1]    : str - Cigar line alignment annotation (M = Align & Seq match, m = Align matcht & Seq mismatch, D = Deletion in ProbeFeature wrt genome, U = Unknown at time of alignment)
  Example    : my $cg = $feature->cigar_string();
  Description: Getter and setter for number of the cigar line attribute for this feature.
  Returntype : str
  Exceptions : None
  Caller     : General
  Status     : High Risk

=cut

sub cigar_string {
  my $self = shift;
  $self->{'cigar_string'} = shift if @_;
  return $self->{'cigar_string'};
}

=head2 probe

  Arg [1]    : Bio::EnsEMBL::Funcgen::Probe - probe
  Example    : my $probe = $feature->probe();
  Description: Getter, setter and lazy loader of probe attribute for
               ProbeFeature objects. Features are retrieved from the database
               without attached probes, so retrieving probe information for a
               feature will involve another query.
  Returntype : Bio::EnsEMBL::Funcgen::Probe
  Exceptions : None
  Caller     : General
  Status     : at risk

=cut

sub probe {
  my $self  = shift;
  my $probe = shift;
  
  if ($probe) {
    if ( !ref $probe || !$probe->isa('Bio::EnsEMBL::Funcgen::Probe') ) {
      throw('Probe must be a Bio::EnsEMBL::Funcgen::Probe object');
    }
    $self->{'probe'} = $probe;
  }

  if ( ! defined $self->{'probe'}){
    $self->{'probe'} = $self->adaptor()->db()->get_ProbeAdaptor()->fetch_by_dbID($self->probe_id());
  }
  return $self->{'probe'};
}

=head2 probe_id

  Example    : my $probe_id = $pfeature->probe_id();
  Description: Getter for the probe db id of the ProbeFeature
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : at risk

=cut

sub probe_id {
  my $self = shift;
  return $self->{'probe_id'} || $self->probe->dbID();
}

=head2 summary_as_hash

  Example       : $probe_feature_summary = $probe_feature->summary_as_hash;
  Description   : Retrieves a textual summary of this ProbeFeature.
  Returns       : Hashref of descriptive strings
  Status        : Intended for internal use (REST)

=cut

sub summary_as_hash {
  my $self = shift;

  my $probe     = $self->probe;
  my $probe_set = defined $probe->probe_set ? $probe->probe_set->name : '';

  return {
    seq_region_name   => $self->seq_region_name,
    start             => $self->seq_region_start,
    end               => $self->seq_region_end,
    strand            => $self->strand,
    feature_type      => 'probe_feature',
    microarray        => $probe->array_chip->name,
    probe_name        => $probe->name,
    probe_set         => $probe_set,
    probe_length      => $probe->length,
   };
}

1;
