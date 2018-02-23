# Ensembl module for Bio::EnsEMBL::Funcgen::BindingMatrix::Converter

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

Bio::EnsEMBL::Funcgen::BindingMatrix::Converter

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::BindingMatrix
Bio::EnsEMBL::Funcgen::DBSQL::BindingMatrixAdaptor
Bio::EnsEMBL::Funcgen::MotifFeature

=cut

package Bio::EnsEMBL::Funcgen::BindingMatrix::Converter;

use strict;
use warnings;
use autodie;
use feature qw(say);

use Bio::EnsEMBL::Utils::Scalar qw( assert_ref check_ref );
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );

use Bio::EnsEMBL::Funcgen::BindingMatrix;
use Bio::EnsEMBL::Funcgen::BindingMatrix::Constants qw ( :all );

sub new {
    my $caller    = shift;
    my $obj_class = ref($caller) || $caller;
    my $self      = {};

    return bless $self, $obj_class;
}


sub from_frequencies_to_probabilities {
    my ( $self, $binding_matrix, $pseudocount ) = @_;

    assert_ref( $binding_matrix, 'Bio::EnsEMBL::Funcgen::BindingMatrix',
        'BindingMatrix' );

    if ( $binding_matrix->unit ne FREQUENCIES ) {
        throw(  'Please supply a binding matrix with '
              . FREQUENCIES
              . ' units instead of '
              . $binding_matrix->unit() );
    }

    my $default_pseudocount = 0.1;
    $pseudocount //= $default_pseudocount;

    my $probabilities = {};

    for (
        my $position = 1 ;
        $position <= $binding_matrix->length() ;
        $position++
      )
    {
        $probabilities->{$position} //= {};
        my $frequency_sum_by_position =
          $self->_get_frequency_sum_by_position( $binding_matrix, $position );

        for my $nucleotide ( @{ $self->_nucleotides() } ) {

            my $frequency =
              $binding_matrix->get_element_by_position_nucleotide( $position,
                $nucleotide );

            $probabilities->{$position}->{$nucleotide} =
              ( $frequency + $pseudocount ) /
              ( $frequency_sum_by_position + 4 * $pseudocount );
        }
    }

    my $probabilities_binding_matrix =
      $self->_convert_BindingMatrix( $binding_matrix, $probabilities,
        PROBABILITIES );

    return $probabilities_binding_matrix;
}

sub from_probabilities_to_weights {
    my $self = shift;
    my $binding_matrix;
    my %expected_frequency;

    (
        $binding_matrix,          $expected_frequency{'A'},
        $expected_frequency{'C'}, $expected_frequency{'G'},
        $expected_frequency{'T'}
      )
      = rearrange(
        [
            'BINDING_MATRIX',       'EXPECTED_FREQUENCY_A',
            'EXPECTED_FREQUENCY_C', 'EXPECTED_FREQUENCY_G',
            'EXPECTED_FREQUENCY_T'
        ],
        @_
      );

    assert_ref( $binding_matrix, 'Bio::EnsEMBL::Funcgen::BindingMatrix',
        'BindingMatrix' );

    if ( $binding_matrix->unit ne PROBABILITIES ) {
        throw(  'Please supply a binding matrix with '
              . PROBABILITIES
              . ' units instead of '
              . $binding_matrix->unit() );
    }

    my $default_expected_frequency = 0.25;

    for my $nucleotide ( keys %expected_frequency ) {
        $expected_frequency{$nucleotide} //= $default_expected_frequency;
    }

    $self->_check_expected_frequencies_are_valid(\%expected_frequency);

    my $weights = {};

    for (
        my $position = 1 ;
        $position <= $binding_matrix->length() ;
        $position++
      )
    {
        $weights->{$position} //= {};

        for my $nucleotide ( @{ $self->_nucleotides() } ) {
            my $probability =
              $binding_matrix->get_element_by_position_nucleotide( $position,
                $nucleotide );

            $weights->{$position}->{$nucleotide} =
              $self->_log2( $probability / $expected_frequency{$nucleotide} );
        }
    }

    my $weights_binding_matrix =
      $self->_convert_BindingMatrix( $binding_matrix, $weights, WEIGHTS );


    return $weights_binding_matrix;
}


sub from_probabilities_to_bits {
    my ( $self, $binding_matrix ) = @_;

    assert_ref( $binding_matrix, 'Bio::EnsEMBL::Funcgen::BindingMatrix',
        'BindingMatrix' );

    if ( $binding_matrix->unit ne PROBABILITIES ) {
        throw(    'Please supply a binding matrix with '
                . PROBABILITIES
                . ' units instead of '
                . $binding_matrix->unit() );
    }

    my $bits = {};

    for (
        my $position = 1;
        $position <= $binding_matrix->length();
        $position++
        )
    {
        $bits->{$position} //= {};
        my $h = $self->_get_h( $binding_matrix, $position );
        my $IC = 2 - $h;

        for my $nucleotide ( @{ $self->_nucleotides() } ) {
            my $probability
                = $binding_matrix->get_element_by_position_nucleotide(
                $position, $nucleotide );

            $bits->{$position}->{$nucleotide} = $probability * $IC;
        }
    }

    my $bits_binding_matrix
        = $self->_convert_BindingMatrix( $binding_matrix, $bits, BITS );

    return $bits_binding_matrix;
}

sub from_frequencies_to_bits {
    my ( $self, $binding_matrix, $pseudocount ) = @_;

    my $probabilities_binding_matrix
        = $self->from_frequencies_to_probabilities($binding_matrix, $pseudocount);

    my $bits_binding_matrix
        = $self->from_probabilities_to_bits($probabilities_binding_matrix);

    return $bits_binding_matrix;
}

sub _get_frequency_sum_by_position {
    my ( $self, $binding_matrix, $position ) = @_;

    assert_ref( $binding_matrix, 'Bio::EnsEMBL::Funcgen::BindingMatrix',
        'BindingMatrix' );

    throw('Must supply a position parameter') if !defined $position;

    my $frequency_sum;

    for my $nucleotide ( @{ $self->_nucleotides() } ) {
        my $frequency
            = $binding_matrix->get_element_by_position_nucleotide( $position,
            $nucleotide );
        $frequency_sum += $frequency;
    }

    return $frequency_sum;
}

sub _log2 {
    my ( $self, $n ) = @_;

    throw('Must supply a parameter') if !defined $n;
    return log($n) / log(2);
}

sub _get_h {
    my ( $self, $binding_matrix, $position ) = @_;

    assert_ref( $binding_matrix, 'Bio::EnsEMBL::Funcgen::BindingMatrix',
        'BindingMatrix' );

    throw('Must specify a position') if !defined $position;

    my $h = 0;

    for my $nucleotide ( @{ $self->_nucleotides() } ) {
        
        my $probability
            = $binding_matrix->get_element_by_position_nucleotide( $position,
            $nucleotide );

        $h -= $probability * $self->_log2($probability);
    }

    return $h;
}


sub _nucleotides { return [ 'A', 'C', 'G', 'T' ]; }

sub _convert_BindingMatrix {
    my ( $self, $binding_matrix, $elements, $unit ) = @_;

    assert_ref( $binding_matrix, 'Bio::EnsEMBL::Funcgen::BindingMatrix',
        'BindingMatrix' );
    throw('Must supply an -elements parameter') if !defined $elements;
    throw('Must supply a -unit parameter')      if !defined $unit;

    return Bio::EnsEMBL::Funcgen::BindingMatrix->new(
        -NAME      => $binding_matrix->name(),
        -SOURCE    => $binding_matrix->source(),
        -THRESHOLD => $binding_matrix->threshold(),
        -ELEMENTS  => $elements,
        -UNIT      => $unit,
        -ASSOCIATED_TRANSCRIPTION_FACTOR_COMPLEXES =>
          $binding_matrix->get_all_associated_TranscriptionFactorComplexes()
    );
}

sub _check_expected_frequencies_are_valid {
    my ( $self, $expected_frequency ) = @_;

    if ( $expected_frequency->{A} +
        $expected_frequency->{C} +
        $expected_frequency->{G} +
        $expected_frequency->{T} != 1 )
    {
        throw(
            'Invalid expected frequencies passed. The sum is not equal to 1.');
    }

    for my $ef ( values $expected_frequency ) {
        if ( $ef == 0 ) {
            throw(
               'Invalid expected frequencies passed. No zero (0) values allowed'
            );
        }
    }
}

1;

