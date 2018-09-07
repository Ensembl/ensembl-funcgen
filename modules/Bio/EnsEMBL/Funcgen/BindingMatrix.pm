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

use List::Util qw(min max);
use Bio::EnsEMBL::Utils::Scalar    qw( assert_ref check_ref assert_integer);
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw deprecate );
use Bio::EnsEMBL::Funcgen::Sequencing::MotifTools qw( parse_matrix_line
                                                      reverse_complement_matrix );
use Bio::EnsEMBL::Funcgen::BindingMatrix::Constants qw ( :all );
require Bio::EnsEMBL::Funcgen::BindingMatrix::Converter;

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
  Exceptions : Throws if name or/and type or/and stable_id not defined
  Caller     : General
  Status     : Medium risk

=cut

sub new {
    my $caller    = shift;
    my $obj_class = ref($caller) || $caller;
    my $self      = $obj_class->SUPER::new(@_);

    my ( $name, $source, $threshold, $elements, $unit,
        $associated_transcription_factor_complexes, $stable_id )
        = rearrange(
        [   'NAME',      'SOURCE',
            'THRESHOLD', 'ELEMENTS', 'UNIT',
            'ASSOCIATED_TRANSCRIPTION_FACTOR_COMPLEXES', 'STABLE_ID'
        ],
        @_
        );

    throw('Must supply a -name parameter')   if !defined $name;
    throw('Must supply a -source parameter') if !defined $source;
    throw('Must supply a -stable_id parameter') if !defined $stable_id;

    $self->{name}      = $name;
    $self->{source}    = $source;
    $self->{stable_id} = $stable_id;
    $self->{threshold} = $threshold if defined $threshold;
    $self->{elements}  = $elements if defined $elements;

    if (defined $unit && $self->_unit_is_valid($unit)){
            $self->{unit} = $unit;
    }
    else {$self->{unit} = FREQUENCIES;}

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

    # $self->_elements();

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
  
  Arg [1]    : String (optional) - Source of binding matrix
  Example    : my $source = $matrix->source;
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

=head2 stable_id
  
  Arg [1]    : String (optional) - Stable Identifier
  Example    : my $stable_id = matrix->stable_id;
  Description: Getter/Setter for source attribute
  Returntype : Scalar - string
  Exceptions : None
  Caller     : General
  Status     : At Risk

=cut

sub stable_id {
  my $self = shift;
  $self->{stable_id} = shift if @_;
  return $self->{stable_id};
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

=head2 get_element_by_position_nucleotide
  
  Arg [1]    : Integer - Position in the matrix
  Arg [2]    : String - Nucleotide of interest (A, C, G or T)
  Example    : my $element = $binding_matrix->get_element_by_position_nucleotide(3,'A');
  Description: Get element value for a particular position and nucleotide
  Returntype : Scalar - numeric value
  Exceptions : Throws if position parameter is not specified
               Throws if nucleotide parameter is not specified
               Throws if positions is out of bounds
               Throws if nucleotide is invalid
  Caller     : General
  Status     : At Risk

=cut

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

=head2 get_TranscriptionFactorComplex_names

  Example    : my $tfc_names =
             : $binding_matrix->get_TranscriptionFactorComplex_names;
  Description: Returns names of all TranscriptionFactorComplexes that are
             : associated with this BindingMatrix
  Returntype : Arrayref of strings
  Exceptions : None

=cut

sub get_TranscriptionFactorComplex_names {

    my ($self) = @_;
    my @names;

    my $transcription_factor_complexes =
      $self->get_all_associated_TranscriptionFactorComplexes;
    for my $transcription_factor_complex ( @{$transcription_factor_complexes} )
    {
        push @names, $transcription_factor_complex->display_name();
    }

    return \@names;
}


sub _min_max_sequence_similarity_score {
    my ($self) = @_;
    if (! $self->{min_sequence_similarity_score}){
        my $converter = Bio::EnsEMBL::Funcgen::BindingMatrix::Converter->new();
        my $weights_binding_matrix =
          $converter->from_frequencies_to_weights( $self );

        
        my @elements_by_position;

        for (my $position = 1; $position<= $weights_binding_matrix->length(); $position++){
            @elements_by_position = values %{$weights_binding_matrix->_elements()->{$position}};
            $self->{min_sequence_similarity_score} += min @elements_by_position;
            $self->{max_sequence_similarity_score} += max @elements_by_position;
        }
    }

    return ($self->{min_sequence_similarity_score}, $self->{max_sequence_similarity_score});
}

=head2 sequence_similarity_score
  
  Arg [1]       : String - sequence of interest
  Example       : $seq_sim_score = $binding_matrix->sequence_similarity_score($seq);
  Description   : Calculates the similarity score of a given sequence
  Returns       : Float
  Exceptions    : Throws if sequence parameter is not specified
                  Throws if sequence contains invalid characters
                  Throws if sequence is not the same length as the binding matrix
  Status        : At risk

=cut

sub sequence_similarity_score {
    my ( $self, $sequence ) = @_;
    my $sequence_similarity_score = 0;

    if ( !$sequence ) {
        throw('Sequence parameter not provided');
    }

    $sequence = uc($sequence);
    $sequence =~ s/\s+//g;

    if ( $sequence =~ /[^ACGT]/ ) {
        throw( 'Sequence ' . $sequence . ' contains invalid characters' );
    }

    if ( CORE::length($sequence) != $self->length() ) {
        throw(  'Specified sequence does not match matrix length!' . "\n"
              . 'Binding Matrix length: '
              . $self->length() . "\n"
              . 'Specified sequence length: '
              . CORE::length($sequence) );
    }

    my $converter = Bio::EnsEMBL::Funcgen::BindingMatrix::Converter->new();
    my $weights_binding_matrix = $converter->from_frequencies_to_weights($self);


    my @bases = split //, $sequence;
    my $position = 1;

    for my $base (@bases) {
        $sequence_similarity_score +=
          $weights_binding_matrix->get_element_by_position_nucleotide(
            $position, $base );
        $position++;
    }

    return $sequence_similarity_score;
}

=head2 relative_sequence_similarity_score
  
  Arg [1]       : String - sequence of interest
  Arg [2]       : Boolean (optional) - Linear scale results (default is log scale)
  Example       : $relative_seq_sim_score = 
                    $binding_matrix->relative_sequence_similarity_score($seq);
  Description   : Calculates the similarity score of a given sequence relative to the
                  optimal site for the matrix.
  Returns       : Integer, between 0 and 1
  Status        : At risk

=cut

sub relative_sequence_similarity_score {
    my ( $self, $sequence, $linear ) = @_;

    my ( $min_sequence_similarity_score, $max_sequence_similarity_score ) =
      $self->_min_max_sequence_similarity_score();

    my $relative_sequence_similarity_score;
    if ($linear) {
        $relative_sequence_similarity_score =
          (
            exp( $self->sequence_similarity_score($sequence) ) -
              exp($min_sequence_similarity_score) ) /
          (
            exp($max_sequence_similarity_score) -
              exp($min_sequence_similarity_score) );
    }
    else {
        $relative_sequence_similarity_score =
          ( $self->sequence_similarity_score($sequence) -
              $min_sequence_similarity_score ) /
          ( $max_sequence_similarity_score - $min_sequence_similarity_score );
    }

    return $relative_sequence_similarity_score;
}

=head2 is_position_informative
  Arg [1]    : Integer - position within the matrix
  Arg [2]    : (Optional) Float - Threshold [0-2] for information content [default is 1.5]
  Example    : $binding_matrix->is_position_informative($position);
  Description: Returns true if position information content is over threshold
  Returntype : Boolean
  Exceptions : Throws if position or threshold is out of bounds
  Caller     : General
  Status     : At Risk
=cut

sub is_position_informative {
    my ( $self, $position, $threshold ) = @_;
    my $is_position_informative = 0;

    if (! $position){
        throw('Position parameter not provided!');
    }

    if ( $position < 1 || $position > $self->length ) {
        throw(  'Position parameter should be between 1 and '
              . $self->length
              . '. You provided: '
              . $position );
    }

    my $default_threshold = 1.5;
    if (! defined $threshold){
        $threshold = $default_threshold;
    }
    elsif ( $threshold < 0 || $threshold > 2 ) {
        throw(  'Threshold parameter should be between 0 and 2. '
              . 'You provided: '
              . $threshold );
    }

    my $converter = Bio::EnsEMBL::Funcgen::BindingMatrix::Converter->new();
    my $bits_binding_matrix =
      $converter->from_frequencies_to_bits( $self );
    
    my $position_bit_score = 0;
    for my $bit_score (values %{$bits_binding_matrix->{elements}->{$position}}){
        $position_bit_score += $bit_score;
    }

    if ($position_bit_score > $threshold){
        $is_position_informative = 1;
    }
    
    return $is_position_informative;

}

sub _max_element {
    my ($self) = @_;

    my $max_element;

    # use Data::Printer; p $self->_elements;
    my @nucleotide_order = ( 'A', 'C', 'G', 'T' );

    for my $nucleotide (@nucleotide_order) {
        for ( my $position = 1 ; $position <= $self->length() ; $position++ ) {
            my $element =
              $self->get_element_by_position_nucleotide( $position,
                $nucleotide );
            if ( $element >= $max_element ) {
                $max_element = $element;
            }
        }
    }

    return $max_element;
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
        stable_id                                 => $self->stable_id(),
        associated_transcription_factor_complexes => \@associated_tfc_names,
        max_element                               => $self->_max_element(),
        elements_string                           => $self->get_elements_as_string,
    };
}


1;
