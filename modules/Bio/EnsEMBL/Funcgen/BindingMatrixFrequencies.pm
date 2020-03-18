# Ensembl module for Bio::EnsEMBL::Funcgen::BindingMatrixFrequencies

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

Bio::EnsEMBL::Funcgen::BindingMatrixFrequencies - A module to represent a BindingMatrixFrequencies. 
In EFG this represents the binding affinities of a Transcription Factor to DNA.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::BindingMatrixFrequencies;

# A C G T frequency matrix
my $freqs = "1126 6975 6741 2506 7171 0 11 13 812 867 899 1332
4583 0 99 1117 0 12 0 0 5637 1681 875 4568
801 181 268 3282 0 0 7160 7158 38 2765 4655 391
661 15 63 266 0 7159 0 0 684 1858 742 880";

my $bmfx = Bio::EnsEMBL::Funcgen::BindingMatrixFrequencies->new(
                                                              - binding_matrix => $binding_matrix,
                                                              - position => $position,
                                                              - nucleotide => $nucleotide
                                                              );


=head1 DESCRIPTION

This class represents information about a BindingMatrixFrequencies, containing the position 
and nucleotide. A BindingMatrixFrequencies is always associated to a BindingMatrix

=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::DBSQL::BindingMatrixFrequenciesAdaptor
Bio::EnsEMBL::Funcgen::DBSQL::BindingMatrix
Bio::EnsEMBL::Funcgen::MotifFeature

=cut

package Bio::EnsEMBL::Funcgen::BindingMatrixFrequencies;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Scalar qw( assert_ref check_ref );
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw deprecate);
use Bio::EnsEMBL::Funcgen::Sequencing::MotifTools qw( parse_matrix_line
    reverse_complement_matrix );

use base qw( Bio::EnsEMBL::Funcgen::Storable );

=head2 new

  Arg [-name]       : Scalar - Name of matrix
  Arg [-analysis]   : Bio::EnsEMBL::Analysis - analysis describing how the matrix was obtained
  Arg [-source]     : String (Mandatory) - A string describing the source of the matrix, i.e. SELEX
  Arg [-threshold]  : Scalar (optional) - Numeric minimum relative affinity for binding sites of this matrix
  Arg [-description]: Scalar (optional) - Descriptiom of matrix
  Example    : my $matrix = Bio::EnsEMBL::Funcgen::BindingMatrixFrequencies->new(
                                                              - binding_matrix => $binding_matrix,
                                                              - position => $position,
                                                              - nucleotide => $nucleotide
                                                                );
  Description: Constructor method for BindingMatrixFrequencies class
  Returntype : Bio::EnsEMBL::Funcgen::BindingMatrixFrequencies
  Exceptions : Throws if name or/and type not defined
  Caller     : General
  Status     : Medium risk

=cut

sub new {
    my $caller    = shift;
    my $obj_class = ref($caller) || $caller;
    my $self      = $obj_class->SUPER::new(@_);

    my ( $position, $nucleotide, $frequency, $binding_matrix )
        = rearrange(
        [ 'POSITION', 'NUCLEOTIDE', 'FREQUENCY', 'BINDING_MATRIX' ], @_ );

    throw('Must supply a -position parameter')   if !defined $position;
    throw('Must supply a -nucleotide parameter') if !defined $nucleotide;
    throw('Must supply a -frequency parameter')  if !defined $frequency;

    if ($binding_matrix) {
        assert_ref( $binding_matrix, 'Bio::EnsEMBL::Funcgen::BindingMatrix',
            'BindingMatrix' );
    }

    $self->{position}       = $position;
    $self->{nucleotide}     = $nucleotide;
    $self->{frequency}      = $frequency;
    $self->{binding_matrix} = $binding_matrix;

    return $self;
}

=head2 get_BindingMatrix

  Example    : my $binding_matrix = $bmf->get_BindingMatrix();
  Description: Getter for the BindingMatrix object
  Returntype : Bio::EnsEMBL::Funcgen::BindingMatrix
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub get_BindingMatrix { return shift->{binding_matrix}; }

=head2 position

  Example    : my $position = $bmf->position();
  Description: Getter for the position attribute
  Returntype : Integer
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub position { return shift->{position}; }

=head2 nucleotide

  Example    : my $nucleotide = $bmf->nucleotide();
  Description: Getter for the nucleotide attribute
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub nucleotide { return shift->{nucleotide}; }

=head2 frequency

  Example    : my $frequency = $bmf->frequency();
  Description: Getter for the frequency attribute
  Returntype : Integer
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub frequency { return shift->{frequency}; }

1;
