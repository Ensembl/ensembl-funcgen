# Ensembl module for Bio::EnsEMBL::Funcgen::TranscriptionFactorComplex

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

package Bio::EnsEMBL::Funcgen::TranscriptionFactorComplex;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Scalar qw( assert_ref check_ref );
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );

use base qw( Bio::EnsEMBL::Funcgen::Storable );

=head2 new

  Arg [-PRODUCTION_NAME] : String - the name which is used internally
  Arg [-DISPLAY_NAME]    : String - the name used on the EnsEMBL browser
  Arg [-COMPONENTS]      : Arrayref of
                           Bio::EnsEMBL::Funcgen::TranscriptionFactor objects
                           All Transcription Factors which are members of
                           this complex

  Example    : $tfc = Bio::EnsEMBL::Funcgen::TrancriptionFactorComplex->new(...);
  Description: Creates a new transcription_factor_complex object
  Returntype : Bio::EnsEMBL::Funcgen::TranscriptionFactorComplex
  Exceptions : Throws if the components parameter contains invalid
               Bio::EnsEMBL::Funcgen::TranscriptionFactor objects
  Caller     : General
  Status     : Stable

=cut

sub new {
    my $caller    = shift;
    my $obj_class = ref($caller) || $caller;
    my $self      = $obj_class->SUPER::new(@_);

    my ( $production_name, $display_name, $components )
        = rearrange( [ 'PRODUCTION_NAME', 'DISPLAY_NAME', 'COMPONENTS' ],
        @_ );

    throw('Must supply a -production_name parameter')
        if !defined $production_name;
    throw('Must supply a -display_name parameter') if !defined $display_name;
    throw('Must supply a -components parameter') if !defined $components;

    for my $transcription_factor ( @{$components} ) {
        assert_ref( $transcription_factor,
            'Bio::EnsEMBL::Funcgen::TranscriptionFactor',
            'TranscriptionFactor' );
    }

    $self->{production_name} = $production_name;
    $self->{display_name}    = $display_name;
    $self->{components}      = $components;

    return $self;
}


=head2 production_name

  Example    : my $production_name = $transcription_factor_complex->production_name();
  Description: Getter for the production name
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub production_name { return shift->{production_name}; }

=head2 display_name

  Example    : my $display_name = $transcription_factor_complex->display_name();
  Description: Getter for the display name
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub display_name { return shift->{display_name}; }

=head2 components

  Example    : my $transcription_factors = $transcription_factor_complex->components();
  Description: Getter for the components of the complex
  Returntype : Arrayref of Bio::EnsEMBL::Funcgen::TranscriptionFactor objects
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub components { return shift->{components}; }

1;
