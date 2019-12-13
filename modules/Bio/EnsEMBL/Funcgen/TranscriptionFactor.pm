# Ensembl module for Bio::EnsEMBL::Funcgen::TranscriptionFactor

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


=cut

package Bio::EnsEMBL::Funcgen::TranscriptionFactor;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Scalar qw( assert_ref check_ref );
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );

use base qw( Bio::EnsEMBL::Funcgen::Storable );

=head2 new

  Arg [-name]           : Scalar (Mandatory) - Name of transcription factor
  Arg [-feature_type]   : Bio::EnsEMBL::Funcgen::FeatureType (optional) -
                          corresponding feature type
  Arg [-gene_stable_id] : String (optional) - Numeric minimum relative affinity
                          for binding sites of this matrix

  Example     : my $matrix = Bio::EnsEMBL::Funcgen::TranscriptionFactor->new(
                             -name  => "CTCF");
  Description : Constructor method for TranscriptionFactor class
  Returntype  : Bio::EnsEMBL::Funcgen::TranscriptionFactor
  Exceptions  : Throws if name not defined
                Throws if feature_type is not a
                Bio::EnsEMBL::Funcgen::FeatureType object
  Caller      : General
  Status      : Medium risk

=cut

sub new {
    my $caller    = shift;
    my $obj_class = ref($caller) || $caller;
    my $self      = $obj_class->SUPER::new(@_);

    my ( $name, $feature_type, $gene_stable_id )
        = rearrange( [ 'NAME', 'FEATURE_TYPE', 'GENE_STABLE_ID' ], @_ );

    if ($feature_type) {
        assert_ref( $feature_type, 'Bio::EnsEMBL::Funcgen::FeatureType',
            'FeatureType' );
    }
    throw('Must supply a -name parameter') if !defined $name;

    $self->{name}           = $name;
    $self->{feature_type}   = $feature_type if $feature_type;
    $self->{gene_stable_id} = $gene_stable_id if $gene_stable_id;

    return $self;
}


=head2 name

  Example    : my $name = $transcription_factor->name();
  Description: Getter for the name
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub name { return shift->{name}; }

=head2 get_FeatureType

  Example    : my $feature_type = $transcription_factor->get_FeatureType();
  Description: Getter for the FeatureType
  Returntype : Bio::EnsEMBL::Funcgen::FeatureType
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub get_FeatureType { return shift->{feature_type}; }

=head2 gene_stable_id

  Example    : my $gene_stable_id = $transcription_factor->gene_stable_id();
  Description: Getter for the gene_stable_id attribute
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub gene_stable_id { return shift->{gene_stable_id}; }

1;
