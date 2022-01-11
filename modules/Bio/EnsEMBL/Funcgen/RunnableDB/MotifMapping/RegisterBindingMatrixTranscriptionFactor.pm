=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::RegisterBindingMatrixTranscriptionFactor

=head1 DESCRIPTION



=head1 LICENSE

    Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
    Copyright [2016-2022] EMBL-European Bioinformatics Institute

    Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

         http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software distributed under the License
    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and limitations under the License.

=head1 CONTACT

    ensembl-funcgen@ebi.ac.uk

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::RegisterBindingMatrixTranscriptionFactor;

use strict;
use base ('Bio::EnsEMBL::Hive::Process');

use Bio::EnsEMBL::Funcgen::DBSQL::TranscriptionFactorAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::TranscriptionFactorComplexAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::BindingMatrixAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::BindingMatrixFrequenciesAdaptor;

use Bio::EnsEMBL::Funcgen::TranscriptionFactor;
use Bio::EnsEMBL::Funcgen::TranscriptionFactorComplex;
use Bio::EnsEMBL::Funcgen::BindingMatrix;
use Bio::EnsEMBL::Funcgen::BindingMatrixFrequencies;


sub run {
    my $self = shift;

    my $species   = $self->param('species');
    my $selex_dir = $self->param_required('selex_dir');

    my $list         = $selex_dir . '/' . $species . '.tsv';
    my $matrices_dir = $selex_dir . '/matrices/';
    my $source       = 'SELEX';


    # -------------------
    # connect to funcgen database
    # -------------------

    my $funcgen_adaptor =
        Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');

    # ---------------
    # read motif list
    # ---------------

    open my $fh, '<', $list;
    my $stable_id_cnt = 1;
    my %already_stored_binding_matrices;

    while (my $current_line = readline $fh) {
        chomp $current_line;

        my $is_header = $current_line =~ /^TF1/;

        #skip header
        if ($is_header) {
            next;
        }

        my @fields = split /\t/, $current_line;

        my $tf1_name     = $fields[0];
        my $tf2_name     = $fields[1];
        my $gene_id1     = $fields[2];
        my $gene_id2     = $fields[3];
        my @matrix_names = split /,/, $fields[8]; #representative motifs

        # -------------------------
        # store to funcgen database
        # -------------------------

        my $first_transcription_factor =
            store_transcription_factor($tf1_name, $gene_id1, $funcgen_adaptor);
        my $second_transcription_factor =
            store_transcription_factor($tf2_name, $gene_id2, $funcgen_adaptor);

        my $transcription_factor_complex =
            store_transcription_factor_complex($first_transcription_factor,
                                               $second_transcription_factor,
                                               $funcgen_adaptor);

        for my $matrix_name (@matrix_names) {
            if (exists $already_stored_binding_matrices{$matrix_name}) {
                my $binding_matrix_adaptor =
                    $funcgen_adaptor->get_adaptor('BindingMatrix');
                my $stored_binding_matrix =
                    $binding_matrix_adaptor->fetch_by_name($matrix_name);
                push @{$stored_binding_matrix
                    ->{associated_transcription_factor_complexes}},
                     $transcription_factor_complex;
                $binding_matrix_adaptor
                    ->_store_binding_matrix_transcription_factor_complex($stored_binding_matrix);
                next;
            }
            store_binding_matrix($matrix_name, $stable_id_cnt, $source,
                                 $funcgen_adaptor, $matrices_dir,
                                 $transcription_factor_complex);
            $already_stored_binding_matrices{$matrix_name} = 1;
            $stable_id_cnt++;
        }
    }

    close $fh;
}

sub store_transcription_factor {
    my ($tf_name, $gene_id, $db_adaptor) = @_;

    if ($tf_name eq 'NA') {return;}

    my $transcription_factor_adaptor
        = $db_adaptor->get_adaptor('TranscriptionFactor');

    my $transcription_factor
        = $transcription_factor_adaptor->fetch_by_name($tf_name);

    if (!$transcription_factor) {

        my $featyre_type_adaptor = $db_adaptor->get_adaptor('FeatureType');
        my $feature_type         =
            $featyre_type_adaptor->fetch_by_name($tf_name);

        $transcription_factor
            = Bio::EnsEMBL::Funcgen::TranscriptionFactor->new(
            -NAME           => $tf_name,
            -FEATURE_TYPE   => $feature_type,
            -GENE_STABLE_ID => $gene_id
        );

        $transcription_factor_adaptor->store($transcription_factor);
        # $logger->info('Transcription Factor ' . $tf_name . ' stored' . "\n",
        #               0, 1);
    }
    # else {
    #     $logger->info(
    #         'Transcription Factor '
    #             . $tf_name
    #             . ' found in DB. Skipping...' . "\n",
    #         0, 1
    #     );
    # }

    return $transcription_factor;
}


sub store_transcription_factor_complex {
    my ($first_transcription_factor, $second_transcription_factor, $db_adaptor)
        = @_;

    my $transcription_factor_complex_adaptor
        = $db_adaptor->get_adaptor('TranscriptionFactorComplex');

    my $production_name = $first_transcription_factor->name;
    my $display_name    = $first_transcription_factor->name;

    my @components;
    push @components, $first_transcription_factor;

    if ($second_transcription_factor) {
        $production_name .= '_' . $second_transcription_factor->name;
        $display_name    .= '::' . $second_transcription_factor->name;
        push @components, $second_transcription_factor;
    }

    my $transcription_factor_complex
        = $transcription_factor_complex_adaptor->fetch_by_production_name(
        $production_name);

    if (!$transcription_factor_complex) {

        $transcription_factor_complex
            = Bio::EnsEMBL::Funcgen::TranscriptionFactorComplex->new(
            -PRODUCTION_NAME => $production_name,
            -DISPLAY_NAME    => $display_name,
            -COMPONENTS      => \@components,
        );

        $transcription_factor_complex_adaptor->store(
            $transcription_factor_complex);

        # $logger->info(
        #     'Transcription Factor Complex '
        #         . $display_name
        #         . ' stored' . "\n",
        #     0, 1
        # );

    }
    # else {
    #     $logger->info(
    #         'Transcription Factor Complex '
    #             . $display_name
    #             . ' found in DB. Skipping...' . "\n",
    #         0, 1
    #     );
    # }

    return $transcription_factor_complex;
}


sub store_binding_matrix {
    my ($matrix_name, $stable_id_cnt, $source, $db_adaptor, $motif_dir,
        $transcription_factor_complex)
        = @_;

    my $binding_matrix_adaptor = $db_adaptor->get_adaptor('BindingMatrix');

    my $binding_matrix = $binding_matrix_adaptor->fetch_by_name($matrix_name);

    if (!$binding_matrix) {

        my $frequencies
            = read_matrix_file($matrix_name, $motif_dir);
        if (!$frequencies) {return;}

        $binding_matrix = Bio::EnsEMBL::Funcgen::BindingMatrix->new(
            -NAME                                      => $matrix_name,
            -THRESHOLD                                 => 4.4,
            -SOURCE                                    => $source,
            -ELEMENTS                                  => $frequencies,
            -ASSOCIATED_TRANSCRIPTION_FACTOR_COMPLEXES =>
            [ $transcription_factor_complex ],
            -STABLE_ID                                 =>
            'ENSPFM' . sprintf("%04d", $stable_id_cnt)
        );

        $binding_matrix_adaptor->store($binding_matrix);

        # $logger->info('Binding Matrix ' . $matrix_name . ' stored' . "\n",
        #               0, 1);
    }
    else {
        $binding_matrix->{associated_transcription_factor_complexes} =
            [ $transcription_factor_complex ];
        $binding_matrix_adaptor
            ->_store_binding_matrix_transcription_factor_complex($binding_matrix);
        # $logger->info(
        #     'Binding Matrix '
        #         . $matrix_name
        #         . ' found in DB. Skipping...' . "\n",
        #     0, 1
        # );
    }
}


sub read_matrix_file {
    my ($matrix_name, $motif_dir) = @_;
    my $line_cnt                  = 0;
    my @frequencies_matrix;

    my $filepath = $motif_dir . '/' . $matrix_name . '.pfm';

    if (!-e $filepath) {
        # $logger->info($filepath . ' not found' . "\n", 0, 1);
        return;
    }

    my @nucleotide_order = ('A', 'C', 'G', 'T');

    open my $fh, '<', $filepath;

    while (my $current_line = readline $fh) {

        chomp $current_line;

        my $is_valid_frequencies_line = $current_line =~ /[\d\s]+/;

        if ($is_valid_frequencies_line) {
            my @frs = split /\s+/, $current_line;

            my $col_cnt = 1;

            for my $fr (@frs) {
                my $frequency =
                    Bio::EnsEMBL::Funcgen::BindingMatrixFrequencies->new(
                        -POSITION   => $col_cnt,
                        -NUCLEOTIDE => $nucleotide_order[$line_cnt],
                        -FREQUENCY  => $fr
                    );

                push @frequencies_matrix, $frequency;

                $col_cnt++;
            }

            $line_cnt++;
        }
    }

    close $fh;

    return \@frequencies_matrix;

}

1;
