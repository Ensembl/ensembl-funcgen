#
# Ensembl module for Bio::EnsEMBL::Funcgen::Peak
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

Bio::EnsEMBL::Peak - A module to represent an enriched feature mapping i.e. a peak call.

=head1 SYNOPSIS

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::Peak;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );

use base (
  'Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality',
  'Bio::EnsEMBL::Funcgen::SetFeature', 
);

=head2 score

  Example    : my $score = $feature->score;
  Description: Getter for the score attribute for this feature. 
  Returntype : String (float)
  Exceptions : None
  Caller     : General
  Status     : Stable

=head2 summit

  Arg [1]    : (optional) int - summit postition
  Example    : my $peak_summit = $feature->summit;
  Description: Getter for the summit attribute for this feature. 
  Returntype : int
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub _constructor_parameters {
  return {
    peak_calling_id   => 'peak_calling_id',
    summit            => 'summit',
    score             => 'score',
    slice             => 'slice',
    seq_region_id     => 'seq_region_id',
    seq_region_start  => 'seq_region_start',
    seq_region_end    => 'seq_region_end',
    seq_region_strand => 'seq_region_strand',
    peak_calling      => 'set_PeakCalling',
  };
}

sub _simple_accessors {
  return [
    { method_name => 'peak_calling_id',    hash_key => '_peak_calling_id',   },
    { method_name => 'summit',             hash_key => 'summit',             },
    { method_name => 'score',              hash_key => 'score',              },
    { method_name => 'start',              hash_key => 'start',              },
    { method_name => 'end',                hash_key => 'end',                },
    { method_name => 'seq_region_id',      hash_key => '_seq_region_id',     },
    { method_name => 'seq_region_start',   hash_key => '_seq_region_start',  },
    { method_name => 'seq_region_end',     hash_key => '_seq_region_end',    },
    { method_name => 'seq_region_strand',  hash_key => '_seq_region_strand', },
    { method_name => 'strand',             hash_key => 'strand',             },
    { method_name => 'slice',              hash_key => 'slice',              },
  ]
}

sub _set_methods {
  return [
    {
      method_name   => 'set_PeakCalling',
      expected_type => 'Bio::EnsEMBL::Funcgen::PeakCalling',
      hash_key      => 'peak_calling',
    },
  ]
}

sub _fetch_methods {
  return [
    {
      method_name             => 'fetch_PeakCalling',
      hash_key                => '_peak_calling',
      get_adaptor_method_name => 'get_PeakCallingAdaptor',
      dbID_method             => 'peak_calling_id',
    },
  ]
}

=head2 display_label

  Example    : my $label = $feature->display_label();
  Description: Getter for the display label of this feature.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub display_label {
    my $self = shift;

    #auto generate here if not set in table
    #need to go with one or other, or can we have both, split into diplay_name and display_label?
    
    if(! $self->{'display_label'}  && $self->adaptor){
      $self->{'display_label'} = $self->feature_type->name()." -";
      $self->{'display_label'} .= " ".$self->epigenome->display_label();
      $self->{'display_label'} .= " Enriched Site";
    }
	
    return $self->{'display_label'};
}

=head2 display_id

  Example    : my $label = $feature->display_id;
  Description: Getter for the display_id of this feature. This was created 
               for generating the display id used in big bed files. Converting
               from bed to bigbed causes problems, if 
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Medium Risk

=cut

sub display_id {
    my $self = shift;

    if(! $self->{'display_id'}  && $self->adaptor){
      $self->{'display_id'} = join '_', 
        $self->feature_type->name(),
        $self->epigenome->production_name(),
        "_Enriched_Site";
    }
    return $self->{'display_id'};
}

=head2 get_underlying_structure

  Example    : my @loci = @{ $af->get_underlying_structure() };
  Description: Returns and array of loci consisting of:
                  (start, (motif_feature_start, motif_feature_end)*, end)
  Returntype : ARRAYREF
  Exceptions : None
  Caller     : General
  Status     : At Risk - This is TFBS specific and could move to TranscriptionFactorFeature

=cut

#This should really be precomputed and stored in the DB to avoid the MF attr fetch
#Need to be aware of projecting here, as these will expire if we project after this method is called

sub get_underlying_structure{
  my $self = shift;

  if(! defined $self->{underlying_structure}){
    my @loci = ($self->start);
	
    foreach my $mf(@{$self->get_associated_MotifFeatures}){
      push @loci, ($mf->start, $mf->end);
    }

    push @loci, $self->end;
	
    $self->{underlying_structure} = \@loci;
  }

  return $self->{underlying_structure};
}

=head2 get_associated_MotifFeatures

  Example    : my @assoc_mfs = @{ $af->get_associated_MotifFeatures };
  Description: Returns and array associated MotifFeature i.e. MotifFeatures
               representing a relevanting PWM/BindingMatrix
  Returntype : ARRAYREF
  Exceptions : None
  Caller     : General
  Status     : At Risk - This is TFBS specific and could move to TranscriptionFactorFeature

=cut

sub get_associated_MotifFeatures{
  my $self = shift;

  if(! defined $self->{assoc_motif_features}){
    my $mf_adaptor = $self->adaptor->db->get_MotifFeatureAdaptor;
		#These need reslicing!
		$self->{assoc_motif_features} = $mf_adaptor->fetch_all_by_Peak($self, $self->slice);
  }

  return $self->{assoc_motif_features};
}

sub SO_term {
  my $self = shift;
  return $self->feature_type->so_accession;
}

=head2 summary_as_hash

  Example       : $segf_summary = $annotf->summary_as_hash;
  Description   : Retrieves a textual summary of this Peak.
  Returns       : Hashref of descriptive strings
  Status        : Intended for internal use (REST)

=cut

sub summary_as_hash {
  my $self = shift;
  my $feature_set = $self->feature_set;

  return
    {
      feature_type     => $self->feature_type->name,
      epigenome        => $self->epigenome->name,
      source           => $feature_set->analysis->logic_name,
      seq_region_name  => $self->seq_region_name,
      start            => $self->seq_region_start,
      end              => $self->seq_region_end,
      description      => $feature_set->display_label,
      strand           => $self->strand,
      summit           => $self->summit,
      score            => $self->score,
    };
}
1;

