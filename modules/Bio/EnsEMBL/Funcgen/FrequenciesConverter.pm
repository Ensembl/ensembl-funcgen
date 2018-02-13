# Ensembl module for Bio::EnsEMBL::Funcgen::FrequenciesConverter

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

Bio::EnsEMBL::Funcgen::FrequenciesConverter

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::BindingMatrix
Bio::EnsEMBL::Funcgen::DBSQL::BindingMatrixAdaptor
Bio::EnsEMBL::Funcgen::MotifFeature

=cut

package Bio::EnsEMBL::Funcgen::FrequenciesConverter;

use strict;
use warnings;
use autodie;
use feature qw(say);

use Bio::EnsEMBL::Utils::Scalar qw( assert_ref check_ref );
use Bio::EnsEMBL::Utils::Exception qw( throw );

sub new {
    my $caller    = shift;
    my $obj_class = ref($caller) || $caller;
    my $self      = {};

    return bless $self, $obj_class;
}

sub get_probabilities {
    my ( $self, $binding_matrix, $pseudocount ) = @_;

    assert_ref( $binding_matrix, 'Bio::EnsEMBL::Funcgen::BindingMatrix',
        'BindingMatrix' );

    my $default_pseudocount = 0.1;
    $pseudocount //= $default_pseudocount;
    my $probabilities = {};
    my $frequencies   = $binding_matrix->frequencies();

    for (
        my $position = 1;
        $position <= $binding_matrix->length();
        $position++
        )
    {
        $probabilities->{$position} //= {};
        my $position_sum
            = $self->_get_position_sum( $frequencies->{$position} );

        for my $nucleotide ( @{ $self->_nucleotides() } ) {
            $probabilities->{$position}->{$nucleotide}
                = ( $frequencies->{$position}->{$nucleotide} + $pseudocount )
                / ( $position_sum + 4 * $pseudocount );
        }
    }
    return $probabilities;
}

sub get_weights {
    my ( $self, $binding_matrix, $expected_frequency_AT,
        $expected_frequency_CG )
        = @_;

    my $default_expected_frequency = 0.25;
    $expected_frequency_AT //= $default_expected_frequency;
    $expected_frequency_CG //= $default_expected_frequency;

    if ( 2 * $expected_frequency_AT + 2 * $expected_frequency_CG != 1 ) {
        throw('Invalid expected frequencies');
    }
    my $weights       = {};
    my $probabilities = $self->get_probabilities($binding_matrix);

    for (
        my $position = 1;
        $position <= $binding_matrix->length();
        $position++
        )
    {
        $weights->{$position} //= {};

        for my $nucleotide ( @{ $self->_nucleotides() } ) {

            my $expected_frequency;
            if ( $nucleotide =~ /^[AT]$/ ) {
                $expected_frequency = $expected_frequency_AT;
            }
            elsif ( $nucleotide =~ /^[CG]$/ ) {
                $expected_frequency = $expected_frequency_CG;
            }

            $weights->{$position}->{$nucleotide}
                = $self->_log2( $probabilities->{$position}->{$nucleotide}
                    / $expected_frequency );
        }
    }
    return $weights;
}

sub get_bits {
    my ( $self, $binding_matrix ) = @_;

    my $bits          = {};
    my $probabilities = $self->get_probabilities($binding_matrix);

    for (
        my $position = 1;
        $position <= $binding_matrix->length();
        $position++
        )
    {
        $bits->{$position} //= {};
        my $h = $self->_get_h( $probabilities, $position );
        my $IC = 2 - $h;
        for my $nucleotide ( @{ $self->_nucleotides() } ) {
            $bits->{$position}->{$nucleotide}
                = $probabilities->{$position}->{$nucleotide} * $IC;
        }
    }

    return $bits;
}

# sub _get_frequencies_hashref {
#     my ( $self, $binding_matrix ) = @_;

#     assert_ref( $binding_matrix, 'Bio::EnsEMBL::Funcgen::BindingMatrix',
#         'BindingMatrix' );

#     my $frequencies = {};

#     for (
#         my $position = 1;
#         $position <= $binding_matrix->length();
#         $position++
#         )
#     {
#         $frequencies->{$position} //= {};
#         for my $nucleotide ( @{ $self->_nucleotides() } ) {
#             $frequencies->{$position}->{$nucleotide}
#                 = $binding_matrix->get_frequency_by_position_nucleotide(
#                 $position, $nucleotide );
#         }
#     }

#     return $frequencies;
# }

sub _get_position_sum {
    my ( $self, $position_frequencies ) = @_;

    throw('Must supply a position_frequencies parameter')
        if !defined $position_frequencies;

    my $position_sum;

    for my $nucleotide ( @{ $self->_nucleotides() } ) {
        $position_sum += $position_frequencies->{$nucleotide};
    }

    return $position_sum;
}

sub _log2 {
    my ( $self, $n ) = @_;

    throw('Must supply a parameter') if !defined $n;
    return log($n) / log(2);
}

sub _get_h {
    my ( $self, $probabilities, $position ) = @_;
    throw('Must supply a probabilities hashref') if !defined $probabilities;
    throw('Must specify a position')             if !defined $position;

    my $h = 0;

    for my $nucleotide ( @{ $self->_nucleotides() } ) {
        $h -= $probabilities->{$position}->{$nucleotide}
            * $self->_log2( $probabilities->{$position}->{$nucleotide} );
    }

    return $h;
}

sub _nucleotides { return [ 'A', 'C', 'G', 'T' ]; }

# sub _to_string {
#     my ( $self, $hashref ) = @_;
#     my $string;

#     my $length = scalar keys %{$hashref};

#     for my $nucleotide ( @{ $self->_nucleotides() } ) {
#         for ( my $position = 1; $position <= $length; $position++ ) {

#             # say $hashref->{$position}->{$nucleotide};
#             $string .= $hashref->{$position}->{$nucleotide} . "\t";
#         }
#         $string .= "\n";
#     }

#     return $string;
# }

1;
