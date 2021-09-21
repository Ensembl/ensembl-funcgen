=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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

package BindingMatrixTest;

use strict;
use warnings;

use Test::More;
use Test::Exception;

use Bio::EnsEMBL::Funcgen::BindingMatrix::Constants qw(:all);

use parent qw(Bio::EnsEMBL::Funcgen::Test);

sub new {
    my $class = shift;
    my $self = $class->SUPER::new(@_);

    my $binding_matrix = $self->{fetched}->[0];
    my $TFCs = $binding_matrix->get_all_associated_TranscriptionFactorComplexes;
    my $n = scalar @{$TFCs} + 1;
    $self->num_method_tests('get_all_associated_TranscriptionFactorComplexes', $n);

    return $self;
}

sub parameters :Test(setup) {
    my $self = shift;

    my %mandatory_constructor_parameters = (
        -name      => 'dummy_name',
        -source    => 'dummy_source',
    );

    $self->{mandatory_constructor_parameters} =
        \%mandatory_constructor_parameters;

    my %optional_constructor_parameters = (
        -threshold => 10,
        -unit      => BITS,
        -stable_id => 'ENSPFM0000',
    );

    my %constructor_parameters = (%mandatory_constructor_parameters,
                                  %optional_constructor_parameters);

    $self->{constructor_parameters} = \%constructor_parameters;

    my %parameters = (
        'sequence'               => 'ACCGGATGTTGCACAAC',
        'invalid_sequence'       => 'BCCGGATGTTGCACAAC',
        'position'               => 10,
        'out_of_bounds_position' => 100,
        'float'                  => 1.5,
        'nucleotide'             => 'C',
        'invalid_nucleotide'     => 'B',
        'threshold'              => 1.40
    );

    $self->{parameters} = \%parameters;
}

sub define_expected :Test(setup) {
    my $self = shift;

    my $elements_string =
        "789\t351\t352\t39\t44\t1312\t634\t501\t2\t3\t910\t66\t227\t246\t"
            . "1312\t1312\t98\t\n322\t868\t1312\t21\t19\t21\t193\t119\t1\t5\t"
            . "7\t1312\t107\t677\t397\t2\t509\t\n523\t444\t90\t1312\t1312\t"
            . "36\t31\t676\t4\t177\t1312\t78\t1312\t12\t12\t19\t273\t\n365\t"
            . "91\t26\t14\t8\t71\t454\t16\t1312\t1312\t164\t261\t64\t635\t"
            . "22\t11\t431\t\n";

    my $elements = {
        '1'  => {
            'C' => '322',
            'T' => '365',
            'G' => '523',
            'A' => '789'
        },
        '2'  => {
            'G' => '444',
            'A' => '351',
            'C' => '868',
            'T' => '91'
        },
        '3'  => {
            'T' => '26',
            'C' => '1312',
            'A' => '352',
            'G' => '90'
        },
        '4'  => {
            'G' => '1312',
            'A' => '39',
            'C' => '21',
            'T' => '14'
        },
        '5'  => {
            'G' => '1312',
            'A' => '44',
            'C' => '19',
            'T' => '8'
        },
        '6'  => {
            'G' => '36',
            'A' => '1312',
            'C' => '21',
            'T' => '71'
        },
        '7'  => {
            'G' => '31',
            'A' => '634',
            'C' => '193',
            'T' => '454'
        },
        '8'  => {
            'A' => '501',
            'G' => '676',
            'T' => '16',
            'C' => '119'
        },
        '9'  => {
            'G' => '4',
            'A' => '2',
            'C' => '1',
            'T' => '1312'
        },
        '10' => {
            'T' => '1312',
            'C' => '5',
            'A' => '3',
            'G' => '177'
        },
        '11' => {
            'A' => '910',
            'G' => '1312',
            'T' => '164',
            'C' => '7'
        },
        '12' => {
            'C' => '1312',
            'T' => '261',
            'G' => '78',
            'A' => '66'
        },
        '13' => {
            'T' => '64',
            'C' => '107',
            'A' => '227',
            'G' => '1312'
        },
        '14' => {
            'T' => '635',
            'C' => '677',
            'A' => '246',
            'G' => '12'
        },
        '15' => {
            'C' => '397',
            'T' => '22',
            'G' => '12',
            'A' => '1312'
        },
        '16' => {
            'A' => '1312',
            'G' => '19',
            'T' => '11',
            'C' => '2'
        },
        '17' => {
            'C' => '509',
            'T' => '431',
            'G' => '273',
            'A' => '98'
        }
    };

    $self->{expected} = {
        'single_element'          => 5,
        'elements_string'         => $elements_string,
        'TFC_name'                => 'ETV2::CEBPD',
        'number_of_TFCs'          => 7,
        'sim_score'               => 21.03,
        'rel_sim_score'           => 0.97,
        'is_position_informative' => 0,
        'unit'                    => FREQUENCIES
    };

    my %summary = (
        name                                      =>
        'ETV2_CEBPD_AAA_TGTAGA40NTGT_RSCGGANNTTGCGYAAN_m1_c2',
        source                                    => 'SELEX',
        threshold                                 => 4.4,
        length                                    => 17,
        elements                                  => $elements,
        unit                                      => $self->{expected}->{unit},
        stable_id                                 => 'ENSPFM0134',
        associated_transcription_factor_complexes => [ 'ETV2::CEBPD',
                                                       'ETV2::TEF',
                                                       'ERF::CEBPD',
                                                       'ELK1::TEF',
                                                       'FLI1::CEBPB',
                                                       'FLI1::CEBPD',
                                                       'ETV5::CEBPD' ],
        max_position_sum                          => 2393,
        elements_string                           => $elements_string
    );

    $self->{expected}->{summary} = \%summary;
}

sub dbIDs_to_fetch {return [134];}

sub getters_setters {
    return [ 'name', 'threshold', 'source', 'stable_id' ];
}

sub constructor_extra :Test(1) {
    my $self = shift;

    my %incomplete_parameters = %{$self->{constructor_parameters}};
    delete $incomplete_parameters{-unit};
    my $new_binding_matrix =
        Bio::EnsEMBL::Funcgen::BindingMatrix->new(%incomplete_parameters);
    is($new_binding_matrix->unit,
       $self->{expected}->{unit},
       '... and Frequencies is the default unit');
}

sub get_element_by_position_nucleotide :Test(6) {
    my $self = shift;

    my $binding_matrix = $self->{fetched}->[0];
    my $pos            = $self->{parameters}->{position};
    my $inv_pos        = $self->{parameters}->{out_of_bounds_position};
    my $float          = $self->{parameters}->{float};
    my $nuc            = $self->{parameters}->{nucleotide};
    my $inv_nuc        = $self->{parameters}->{invalid_nucleotide};

    is($binding_matrix->get_element_by_position_nucleotide($pos, $nuc),
       $self->{expected}->{single_element},
       'get_element_by_position_nucleotide() works'
    );

    throws_ok {
        $binding_matrix->get_element_by_position_nucleotide(undef, $nuc);
    }
        qr/Must supply a position parameter/,
        "... and exception is thrown when position parameter is missing";

    throws_ok {
        $binding_matrix->get_element_by_position_nucleotide($pos, undef);
    }
        qr/Must supply a nucleotide parameter/,
        "... and exception is thrown when nucleotide parameter is missing";

    throws_ok {
        $binding_matrix->get_element_by_position_nucleotide($inv_pos, $nuc);
    }
        qr/The -position parameter has to be an integer between 1 and/,
        "... and exception is thrown when position parameter is out of bounds";

    throws_ok {
        $binding_matrix->get_element_by_position_nucleotide($float, $nuc);
    }
        qr/Attribute position was a number but not an Integer/,
        "... and exception is thrown when position parameter is not an integer";

    throws_ok {
        $binding_matrix->get_element_by_position_nucleotide($pos, $inv_nuc);
    }
        qr/Supplied nucleotide .* not valid/,
        "... and exception is thrown when nucleotide parameter is not valid";
}

sub get_elements_as_string :Test(1) {
    my $self = shift;

    my $binding_matrix = $self->{fetched}->[0];

    is($binding_matrix->get_elements_as_string,
       $self->{expected}->{elements_string},
       'get_elements_as_string() works'
    );
}

sub get_all_associated_TranscriptionFactorComplexes :Test(no_plan) {
    my $self = shift;

    my $binding_matrix = $self->{fetched}->[0];

    my $TFCs = $binding_matrix->get_all_associated_TranscriptionFactorComplexes;

    for my $TFC (@{$TFCs}) {
        isa_ok($TFC, 'Bio::EnsEMBL::Funcgen::TranscriptionFactorComplex');
    }

    is(scalar @{$TFCs},
       $self->{expected}->{number_of_TFCs},
       'get_all_associated_TranscriptionFactorComplexes() works');
}

sub get_TranscriptionFactorComplex_names :Test(1) {
    my $self = shift;

    my $binding_matrix = $self->{fetched}->[0];

    my $TFC_names = $binding_matrix->get_TranscriptionFactorComplex_names;

    is($TFC_names->[0],
       $self->{expected}->{TFC_name},
       'get_TranscriptionFactorComplex_names() works');
}

sub sequence_similarity_score :Test(4) {
    my $self = shift;

    my $binding_matrix = $self->{fetched}->[0];
    my $sequence       = $self->{parameters}->{sequence};
    my $score          = $binding_matrix->sequence_similarity_score($sequence);
    $score             = sprintf("%.2f", $score);

    is($score,
       $self->{expected}->{sim_score},
       'sequence_similarity_score() works');

    throws_ok {
        $binding_matrix->sequence_similarity_score();
    }
        qr/Sequence parameter not provided/,
        "... and exception is thrown when sequence parameter is missing";

    my $invalid_sequence = $self->{parameters}->{invalid_sequence};
    throws_ok {
        $binding_matrix->sequence_similarity_score($invalid_sequence);
    }
        qr/contains invalid characters/,
        "... and exception is thrown when sequence contains invalid characters";

    my $shortened_sequence = chop $sequence;
    throws_ok {
        $binding_matrix->sequence_similarity_score($shortened_sequence);
    }
        qr/Specified sequence does not match matrix length/,
        "... and exception is thrown when sequence length is incorrect";
}

sub relative_sequence_similarity_score :Test(1) {
    my $self = shift;

    my $binding_matrix = $self->{fetched}->[0];
    my $sequence       = $self->{parameters}->{sequence};
    my $score          =
        $binding_matrix->relative_sequence_similarity_score($sequence);
    $score = sprintf("%.2f", $score);

    is($score,
       $self->{expected}->{rel_sim_score},
       'relative_sequence_similarity_score() works');

    #TODO test with optional linear parameter set to 1
}

sub is_position_informative :Test(4) {
    my $self = shift;

    my $binding_matrix = $self->{fetched}->[0];
    my $position       = $self->{parameters}->{position};
    my $threshold      = $self->{parameters}->{threshold};

    is($binding_matrix->is_position_informative($position),
       $self->{expected}->{is_position_informative},
       'is_position_informative() works with default threshold');

    is($binding_matrix->is_position_informative($position, $threshold),
       !$self->{expected}->{is_position_informative},
       '... and works with a user defined threshold too');

    throws_ok {
        $binding_matrix->is_position_informative();
    }
        qr/Position parameter not provided/,
        "... and exception is thrown when position parameter is missing";

    my $invalid_position = $self->{parameters}->{out_of_bounds_position};
    throws_ok {
        $binding_matrix->is_position_informative($invalid_position);
    }
        qr/Position parameter should be between/,
        "... and exception is thrown when position parameter is out of bounds";

    #TODO test threshold related exception throws
}

1;