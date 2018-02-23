# Ensembl module for Bio::EnsEMBL::Funcgen::BindingMatrix

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

Bio::EnsEMBL::Funcgen::BindingMatrix - A module to represent a BindingMatrix.
In EFG this represents the binding affinities of a Transcription Factor to DNA.

=head1 SYNOPSIS

=head1 DESCRIPTION

This class represents information about a BindingMatrix

=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::DBSQL::BindingMatrixAdaptor
Bio::EnsEMBL::Funcgen::MotifFeature

=cut

package Bio::EnsEMBL::Funcgen::BindingMatrix;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Scalar    qw( assert_ref check_ref assert_integer);
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Funcgen::Sequencing::MotifTools qw( parse_matrix_line
                                                      reverse_complement_matrix );
use Bio::EnsEMBL::Funcgen::BindingMatrix::Constants qw ( :all );

use base qw( Bio::EnsEMBL::Funcgen::Storable );

=head2 new

  Arg [-name]       : Scalar - Name of matrix
  Arg [-analysis]   : Bio::EnsEMBL::Analysis - analysis describing how the matrix was obtained
  Arg [-source]     : String (Mandatory) - A string describing the source of the matrix, i.e. SELEX
  Arg [-threshold]  : Scalar (optional) - Numeric minimum relative affinity for binding sites of this matrix
  Arg [-description]: Scalar (optional) - Descriptiom of matrix
  Example    : my $matrix = Bio::EnsEMBL::Funcgen::BindingMatrix->new(
                                                               -name  => "MA0122.1",
                                                               -analysis => $analysis,
                                                               -description => "Jaspar Matrix",
                                                                );
  Description: Constructor method for BindingMatrix class
  Returntype : Bio::EnsEMBL::Funcgen::BindingMatrix
  Exceptions : Throws if name or/and type not defined
  Caller     : General
  Status     : Medium risk

=cut

sub new {
    my $caller    = shift;
    my $obj_class = ref($caller) || $caller;
    my $self      = $obj_class->SUPER::new(@_);

    my ( $name, $source, $threshold, $elements, $unit,
        $associated_transcription_factor_complexes )
        = rearrange(
        [   'NAME',      'SOURCE',
            'THRESHOLD', 'ELEMENTS', 'UNIT',
            'ASSOCIATED_TRANSCRIPTION_FACTOR_COMPLEXES'
        ],
        @_
        );

    throw('Must supply a -name parameter')   if !defined $name;
    throw('Must supply a -source parameter') if !defined $source;

    $self->{name}      = $name;
    $self->{source}    = $source;
    $self->{threshold} = $threshold if defined $threshold;
    $self->{elements}  = $elements if defined $elements;

    if (defined $unit && $self->_unit_is_valid($unit)){
            $self->{unit} = $unit;
    }

    if ( defined $associated_transcription_factor_complexes ) {
        for my $complex ( @{$associated_transcription_factor_complexes} ) {
            assert_ref( $complex,
                'Bio::EnsEMBL::Funcgen::TranscriptionFactorComplex',
                'TranscriptionFactorComplex' );
        }

        $self->{associated_transcription_factor_complexes} =
          $associated_transcription_factor_complexes;
    }

    return $self;
}

=head2 name

  Example    : my $name = $matrix->name();
  Description: Getter for the name attribute
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub name { return shift->{name}; }

=head2 unit

  Example    : my $unit = $matrix->unit();
  Description: Getter/Setter for the unit attribute
  Returntype : String
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub unit {
    my ( $self, $unit ) = @_;

    $self->_elements();

    if ($unit && $self->_unit_is_valid($unit)) {
        $self->{unit} = $unit;
    }

    return $self->{unit};
}

=head2 threshold

  Arg [1]    : Scalar (optional) - Numeric threshold
  Example    : if($score >= $matrix->threshold) { # Do something here }
  Description: Getter/setter for threshold attribute
  Returntype : Scalar - numeric
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub threshold {
  my $self = shift;
  $self->{threshold} = shift if @_;
  return $self->{threshold};
}

=head2 source

  Example    : print '>'.$matrix->name."\n".$matrix->source;
  Description: Getter/Setter for source attribute
  Returntype : Scalar - string
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub source {
  my $self = shift;
  $self->{source} = shift if @_;
  return $self->{source};
}

sub _unit_is_valid {
    my ( $self, $unit ) = @_;

    my $valid_units = VALID_UNITS;

    if ( grep $_ eq $unit, @{$valid_units} ) {
        return 1;
    }
    else {
        throw(    $unit
                . ' is not a valid BindingMatrix unit. List of valid units : '
                . join( ",", @{$valid_units} ) );
    }
}

sub _elements {
    my ($self) = @_;

    if ( !$self->{elements} ) {

        my $binding_matrix_frequencies_adaptor
            = $self->adaptor->db()->get_adaptor('BindingMatrixFrequencies');

        my $binding_matrix_frequencies
            = $binding_matrix_frequencies_adaptor->fetch_all_by_BindingMatrix(
            $self);

        for my $bmf ( @{$binding_matrix_frequencies} ) {
            $self->{elements}->{ $bmf->position() }->{ $bmf->nucleotide() }
                = $bmf->frequency();
        }

        $self->unit(FREQUENCIES);
    }

    return $self->{elements};
}



sub get_element_by_position_nucleotide {
    my ( $self, $position, $nucleotide ) = @_;

    throw('Must supply a position parameter')   if !defined $position;
    throw('Must supply a nucleotide parameter') if !defined $nucleotide;

    my %valid_nucleotides = ( 'A' => 1, 'C' => 1, 'G' => 1, 'T' => 1 );

    if (!(      assert_integer( $position, 'position' )
            and 1 <= $position
            and $position <= $self->length
        )
        )
    {
        throw( 'The -position parameter has to be an integer between 1 and '
                . $self->length() );
    }

    if ( !$valid_nucleotides{$nucleotide} ) {
        throw('Supplied nucleotide not valid');
    }

    return $self->_elements()->{$position}->{$nucleotide};
}

sub get_elements_as_string {
    my ($self) = @_;

    my $elements_string;

    my @nucleotide_order = ( 'A', 'C', 'G', 'T' );

    for my $nucleotide (@nucleotide_order) {
        for ( my $position = 1; $position <= $self->length(); $position++ ) {
            my $element
                = $self->get_element_by_position_nucleotide( $position,
                $nucleotide );
            $elements_string .= $element . "\t";
        }
        $elements_string .= "\n";
    }

    return $elements_string;
}


=head2 length

  Example    : my $length = $bm->length;
  Description: Returns the length of the the matrix
  Returntype : Scalar - numeric (int)
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub length {
    my $self = shift;

    if ( !$self->{length} ) {
        $self->{length} = scalar keys %{ $self->_elements() };
    }

    return $self->{length};
}

=head2 get_all_associated_TranscriptionFactorComplexes

  Example    : my $associated_tfcs =
             : $binding_matrix->get_all_associated_TranscriptionFactorComplexes;
  Description: Returns all TranscriptionFactorComplexes that are associated with
             : this BindingMatrix
  Returntype : Arrayref of Bio::EnsEMBL::Funcgen::TranscriptionFactorComplex objects
  Exceptions : None

=cut

sub get_all_associated_TranscriptionFactorComplexes {
    my ($self) = @_;

    if ( !$self->{associated_transcription_factor_complexes} ) {
        my $transcription_factor_complex_adaptor =
          $self->adaptor->db->get_adaptor('TranscriptionFactorComplex');

        $self->{associated_transcription_factor_complexes} =
          $transcription_factor_complex_adaptor->fetch_all_by_BindingMatrix(
            $self);
    }

    return $self->{associated_transcription_factor_complexes};
}

=head2 summary_as_hash

  Example       : $binding_matrix_summary = $binding_matrix->summary_as_hash;
  Description   : Retrieves a textual summary of this BindingMatrix.
  Returns       : Hashref of descriptive strings
  Status        : Intended for internal use (REST)

=cut

sub summary_as_hash {
    my $self = shift;

    my @associated_tfc_names;

    for my $complex (
        @{ $self->get_all_associated_TranscriptionFactorComplexes() } )
    {
        push @associated_tfc_names, $complex->display_name();
    }

    return {
        name                                      => $self->name(),
        source                                    => $self->source(),
        threshold                                 => $self->threshold(),
        length                                    => $self->length(),
        elements                                  => $self->_elements(),
        unit                                      => $self->unit(),
        associated_transcription_factor_complexes => \@associated_tfc_names

    };
}


1;
