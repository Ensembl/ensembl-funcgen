#
# Ensembl module for Bio::EnsEMBL::Funcgen::MirnaTargetFeature
#

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::MirnaTargetFeature - A module to represent an externally curated feature
mapping from an external_db.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::MirnaTargetFeature;

my $feature = Bio::EnsEMBL::Funcgen::MirnaTargetFeature->new(
                -SLICE          => $chr_1_slice,
                -START          => 1_000_000,
                -END            => 1_000_024,
                -STRAND         => -1,
                -Analysis       => Bio::EnsEMBL::Analysis,
                -FEATURE_TYPE   => Bio::EnsEMBL::Funcgen::FeatureType,
                -GENE_STABLE_ID => ENSG00000139618
                -DISPLAY_LABEL  => $text,
                                   );




=head1 DESCRIPTION

An MirnaTargetFeature object represents the genomic placement of an externally curated
feature from and DB external to Ensembl.

=cut

package Bio::EnsEMBL::Funcgen::MirnaTargetFeature;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );

# use base qw(Bio::EnsEMBL::Funcgen::SetFeature);
use base qw(Bio::EnsEMBL::Feature Bio::EnsEMBL::Funcgen::Storable);


=head2 new

  Arg [-FEATURE_TYPE]  : Bio::EnsEMBL::Funcgen::FeatureType
  Arg [-ANALYSIS]      : Bio::EnsEMBL::Analysis
  Arg [-SLICE]         : Bio::EnsEMBL::Slice - The slice on which this feature is.
  Arg [-START]         : int - The start coordinate of this feature relative to the start of the slice
		                          it is sitting on. Coordinates start at 1 and are inclusive.
  Arg [-END]           : int - The end coordinate of this feature relative to the start of the slice
	                            it is sitting on. Coordinates start at 1 and are inclusive.
  Arg [-STRAND]        : int - The orientation of this feature. Valid values are 1, -1 and 0.
  Arg [-DISPLAY_LABEL] : string - Display label for this feature
  Arg [-GENE_STABLE_ID]: string - Ensembl Gene Stable ID (ENSG)
  Arg [-dbID]          : (optional) int - Internal database ID.
  Arg [-ADAPTOR]       : (optional) Bio::EnsEMBL::DBSQL::BaseAdaptor - Database adaptor.
  Example              : my $feature = Bio::EnsEMBL::Funcgen::MirnaTargetFeature->new(
                            -SLICE          => $chr_1_slice,
                            -START          => 1_000_000,
                            -END            => 1_000_024,
                            -STRAND         => -1,
                            -Analysis       => Bio::EnsEMBL::Analysis,
                            -FEATURE_TYPE   => Bio::EnsEMBL::Funcgen::FeatureType,
                            -GENE_STABLE_ID => ENSG00000139618
                            -DISPLAY_LABEL  => $text,

                                               );


  Description: Constructor for MirnaTargetFeature objects.
  Returntype : Bio::EnsEMBL::Funcgen::MirnaTargetFeature
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub new {
  my $caller = shift;

  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new(@_);

    my ($feature_type, $analysis, $gene_stable_id, $accession, $evidence,
        $method, $supporting_information, $display_label) =
        rearrange([ 'FEATURE_TYPE', 'ANALYSIS', 'GENE_STABLE_ID', 'ACCESSION',
                    'EVIDENCE', 'METHOD', 'SUPPORTING_INFORMATION',
                    'DISPLAY_LABEL' ], @_);

    for my $var ($feature_type, $analysis, $gene_stable_id, $accession,
                 $evidence, $method, $supporting_information, $display_label) {
        throw 'Must supply a mandatory parameter' unless defined($var) and length $var;
    }

  $self->{feature_type}           = $feature_type;
  $self->{analysis}               = $analysis;
  $self->{gene_stable_id}         = $gene_stable_id;
  $self->{accession}              = $accession;
  $self->{evidence}               = $evidence;
  $self->{method}                 = $method;
  $self->{supporting_information} = $supporting_information;
  $self->{display_label}          = $display_label;

  return $self;
}

sub get_FeatureType {
  return $_[0]->{'feature_type'};
}


=head2 display_label

  Example    : my $label = $feature->display_label();
  Description: Getter for the display label of this feature.
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Medium risk

=cut

sub display_label {
  my $self = shift;

  if(! $self->{'display_label'}  && $self->adaptor){
	 $self->{'display_label'}  = $self->feature_type->name();
  }

  return $self->{'display_label'};
}

=head2 accession

  Arg [1]    : (optional) int - stable_id e.g 1
  Example    : my $acc = $mirna_target_feature->accession();
  Description: Getter for the accession attribute for this MirnaTargetFeature.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub accession {
  return $_[0]->{'accession'};
}

=head2 evidence

  Arg [1]    : (optional) int - stable_id e.g 1
  Example    : my $acc = $mirna_target_feature->evidence();
  Description: Getter for the evidence attribute for this MirnaTargetFeature.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub evidence {
  return $_[0]->{'evidence'};
}

=head2 method

  Arg [1]    : (optional) int - stable_id e.g 1
  Example    : my $acc = $mirna_target_feature->method();
  Description: Getter for the method attribute for this MirnaTargetFeature.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub method {
  return $_[0]->{'method'};
}

=head2 gene_stable_id

  Example    : my $gene_stable_idID = $mirna_target_feature->gene_stable_id();
  Description: Getter for the linked Gene Ensembl StableID for this MirnaTargetFeature.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub gene_stable_id {
  return $_[0]->{'gene_stable_id'};
}


=head2 supporting_information

  Arg [1]    : (optional) int - stable_id e.g 1
  Example    : my $acc = $mirna_target_feature->supporting_information();
  Description: Getter/Setter for the supporting_information attribute for this MirnaTargetFeature.
  Returntype : string
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub supporting_information {
  my ($self, $info) = @_;

  if(defined $info){
    $self->{supporting_information} = $info;
  
  }
  return $self->{'supporting_information'};
}


1;

