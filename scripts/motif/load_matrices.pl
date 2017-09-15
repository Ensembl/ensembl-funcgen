#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use autodie;
use feature qw(say);

use Getopt::Long;
use Data::Printer;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Utils::Logger;

use Bio::EnsEMBL::Funcgen::DBSQL::TranscriptionFactorAdaptor;
# use Bio::EnsEMBL::Funcgen::DBSQL::FeatureTypeAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::TranscriptionFactorComplexAdaptor;
use Bio::EnsEMBL::Funcgen::TranscriptionFactor;
# use Bio::EnsEMBL::Funcgen::FeatureType;
use Bio::EnsEMBL::Funcgen::TranscriptionFactorComplex;

main();

sub main {

    my ( $motif_list, $motif_dir, $source, $overwrite );
    my $logger = Bio::EnsEMBL::Utils::Logger->new();

    GetOptions(
        'motif_list=s' => \$motif_list,
        'motif_dir=s'  => \$motif_dir,
        'source=s'     => \$source,
        'overwrite=i'  => \$overwrite,
    );

    if ( !$motif_list ) {
        $logger->error(
            'Please specify the tab delimited motif list using -motif_list option',
            0, 1
        );
    }

    if ( !-e $motif_list ) {
        $logger->error(
            'Motif list not found. Please make sure that file exists',
            0, 1 );
    }

    if ( !$motif_dir ) {
        $logger->error(
            'Please specify the motif directory using -motif_dir option',
            0, 1 );
    }

    if ( !-d $motif_dir ) {
        $logger->error(
            'Motif directory not found. Please make sure that it exists',
            0, 1 );
    }

    if ( !$source ) {
        $logger->error(
            'Please specify the motif source (i.e. SELEX) using -source option',
            0, 1
        );
    }

    # -------------------
    # connect to funcgen database
    # -------------------

    my $db_adaptor = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(

        -host   => 'mysql-ens-reg-prod-1.ebi.ac.uk',
        -user   => 'ensadmin',
        -pass   => 'ensembl',
        -port   => 4526,
        -dbname => 'ilavidas_motif_homo_sapiens_funcgen_89_38'
    );

    # ---------------
    # read motif list
    # ---------------

    open my $fh, '<', $motif_list;

    while ( readline $fh ) {

        chomp;

        #skip header
        if (/^TF1/) {
            next;
        }

        my @fields = split /\t/;

        my $tf1_name   = $fields[0];
        my $tf2_name   = $fields[1];
        my $gene_id1   = $fields[2];
        my $gene_id2   = $fields[3];
        my $motif_name = $fields[8];    #representative motif

        # -------------------------
        # store to funcgen database
        # -------------------------
        my $first_transcription_factor
            = store_transcription_factor( $tf1_name, $gene_id1, $db_adaptor,
            $logger );
        my $second_transcription_factor
            = store_transcription_factor( $tf2_name, $gene_id2, $db_adaptor,
            $logger );

        store_transcription_factor_complex( $first_transcription_factor,
            $second_transcription_factor, $db_adaptor, $logger );

        # store_binding_matrix( $motif, $db_adaptor );
        # die;
    }

    close $fh;

}

sub store_transcription_factor {
    my ( $tf_name, $gene_id, $db_adaptor, $logger ) = @_;

    if ( $tf_name eq 'NA' ) { return; }

    my $transcription_factor_adaptor
        = $db_adaptor->get_adaptor('TranscriptionFactor');

    my $transcription_factor
        = $transcription_factor_adaptor->fetch_by_name($tf_name);

    if ( !$transcription_factor ) {

        my $featyre_type_adaptor = $db_adaptor->get_adaptor('FeatureType');
        my $feature_type = $featyre_type_adaptor->fetch_by_name($tf_name);

        $transcription_factor
            = Bio::EnsEMBL::Funcgen::TranscriptionFactor->new(
            -NAME           => $tf_name,
            -FEATURE_TYPE   => $feature_type,
            -GENE_STABLE_ID => $gene_id
            );

        $transcription_factor_adaptor->store($transcription_factor);
        $logger->info( 'Transcription Factor ' . $tf_name . ' stored' . "\n",
            0, 1 );
    }
    else {
        $logger->info(
            'Transcription Factor '
                . $tf_name
                . ' found in DB. Skipping...' . "\n",
            0, 1
        );
    }

    return $transcription_factor;
}


sub store_transcription_factor_complex {
    my ( $first_transcription_factor, $second_transcription_factor,
        $db_adaptor, $logger )
        = @_;

    my $transcription_factor_complex_adaptor
        = $db_adaptor->get_adaptor('TranscriptionFactorComplex');

    my $production_name = $first_transcription_factor->name;
    my $display_name = $first_transcription_factor->name;

    if ($second_transcription_factor) {
        $production_name .= '_' . $second_transcription_factor->name;
        $display_name    .= '::' . $second_transcription_factor->name;
    }

    my $transcription_factor_complex
        = $transcription_factor_complex_adaptor->fetch_by_production_name(
        $production_name);

    if ( !$transcription_factor_complex ) {

        $transcription_factor_complex
            = Bio::EnsEMBL::Funcgen::TranscriptionFactorComplex->new(
            -PRODUCTION_NAME => $production_name,
            -DISPLAY_NAME    => $display_name,
            );

        $transcription_factor_complex_adaptor->store(
            $transcription_factor_complex);

    }

}



sub store_transcription_factor_complex_composition {
    my ( $transcription_factor_complex, $transcription_factor, $db_adaptor )
        = @_;

    my $transcription_factor_complex_composition_adaptor
        = $db_adaptor->get_adaptor('TranscriptionFactorComplexComposition');

    my $transcription_factor_complex_composition
        = Bio::EnsEMBL::Funcgen::TranscriptionFactorComplexComposition->new(
        -TRANSCRIPTION_FACTOR_COMPLEX => $transcription_factor_complex,
        -TRANSCRIPTION_FACTOR         => $transcription_factor,
        );

    $transcription_factor_complex_composition_adaptor->store(
        $transcription_factor_complex_composition);

}

sub store_binding_matrix{

}

# sub create_matrix_object {
#     my ( $logger, $motif_name, $tf_name, $source, $db_adaptor ) = @_;

#     my $analysis_adaptor = $db_adaptor->get_adaptor('analysis');
#     my $analysis         = $analysis_adaptor->fetch_by_logic_name('SELEX');

#     my $feature_type_adaptor = $db_adaptor->get_adaptor('featuretype');
#     my $feature_type         = $feature_type_adaptor->fetch_by_name($tf_name);

#     if ( !$feature_type ) {
#         $logger->warning(
#             'Skipping '
#                 . $motif_name
#                 . '. Transcription Factor '
#                 . $tf_name
#                 . ' not found in database.',
#             0, 1
#         );
#         return;
#     }

#     my $binding_matrix = Bio::EnsEMBL::Funcgen::BindingMatrix->new(
#         'name'         => $motif_name,
#         'analysis'     => $analysis,
#         'source'       => $source,
#         'feature_type' => $feature_type,
#     );

#     return $binding_matrix;
# }

# sub create_frequencies_objects {
#     my ( $logger, $frequencies_matrix, $db_adaptor, $binding_matrix ) = @_;

#     my @binding_matrix_frequencies;
#     my $binding_matrix_adaptor = $db_adaptor->get_adaptor('binding_matrix');

#     for my $nucleotide ( keys %{$frequencies_matrix} ) {
#         for (
#             my $position = 0;
#             $position < scalar @{ $frequencies_matrix->{$nucleotide} };
#             $position++
#             )
#         {
#             my $bmf = Bio::EnsEMBL::Funcgen::BindingMatrixFrequencies->new(
#                 'binding_matrix' => $binding_matrix,
#                 'position'       => $position + 1,
#                 'nucleotide'     => $nucleotide,
#                 'frequency' =>
#                     $frequencies_matrix->{$nucleotide}->[$position],
#             );

#             push @binding_matrix_frequencies, $bmf;
#         }
#     }

#     return \@binding_matrix_frequencies;
# }

# opendir my ($dirh), $motif_dir;

# while ( my $motif_file = readdir $dirh ) {

#     if ( $motif_file eq "." || $motif_file eq ".." ) { next; }

#     my ( $motif_name, $tf_name );

#     if ( $motif_file =~ /(.*)\.pfm/ ) {
#         $motif_name = $1;

#         if ( $motif_name =~ /(.*?)_/ ) {
#             $tf_name = $1;
#         }

#         my $binding_matrix
#             = create_matrix_object( $logger, $motif_name, $tf_name,
#             $source, $db_adaptor );

#         if ( !$binding_matrix ) {
#             next;
#         }

#         my $binding_matrix_adaptor
#             = $db_adaptor->get_adaptor('binding_matrix');
#         $binding_matrix_adaptor->store($binding_matrix);

#         open my $fh, '<', $motif_dir . '/' . $motif_file;
#         my $frequencies_matrix;
#         my $line_cnt = 0;

#         while ( readline $fh ) {
#             my @nucleotide_order = ( 'A', 'C', 'G', 'T' );

#             if (/[\d\s]+/) {
#                 my @frs = split /\s+/;
#                 $frequencies_matrix->{ $nucleotide_order[$line_cnt] }
#                     = \@frs;
#                 $line_cnt++;
#             }
#         }

#         close $fh;

#         my $binding_matrix_frequencies = create_frequencies_objects(
#             $logger,     $frequencies_matrix,
#             $db_adaptor, $binding_matrix
#         );

#         my $binding_matrix_frequencies_adaptor = $db_adaptor->get_adaptor(
#             'binding_matrix_frequencies_adaptor');

#         for my $bmf ( @{binding_matrix_frequencies} ) {
#             $binding_matrix_frequencies_adaptor->store($bmf);
#         }

#     }
