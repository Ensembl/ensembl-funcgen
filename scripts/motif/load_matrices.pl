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
use Bio::EnsEMBL::Funcgen::DBSQL::BindingMatrixAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::BindingMatrixFrequenciesAdaptor;
use Bio::EnsEMBL::Funcgen::BindingMatrix;
use Bio::EnsEMBL::Funcgen::BindingMatrixFrequencies;
use Bio::EnsEMBL::Funcgen::DBSQL::FeatureTypeAdaptor;
use Bio::EnsEMBL::Funcgen::FeatureType;

main();

sub main {

    my ( $motif_dir, $source, $overwrite );
    my $logger = Bio::EnsEMBL::Utils::Logger->new();

    GetOptions(
        'motif_dir=s' => \$motif_dir,
        'source=s'    => \$source,
        'overwrite=i' => \$overwrite,
    );

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

    # -------------------
    # read motif directory
    # -------------------

    opendir my ($dirh), $motif_dir;

    while ( my $motif_file = readdir $dirh ) {

        if ( $motif_file eq "." || $motif_file eq ".." ) { next; }

        my ( $motif_name, $tf_name );

        if ( $motif_file =~ /(.*)\.pfm/ ) {
            $motif_name = $1;

            if ( $motif_name =~ /(.*?)_/ ) {
                $tf_name = $1;
            }

            my $binding_matrix
                = create_matrix_object( $logger, $motif_name, $tf_name,
                $source, $db_adaptor );

            if ( !$binding_matrix ) {
                next;
            }

            my $binding_matrix_adaptor
                = $db_adaptor->get_adaptor('binding_matrix');
            $binding_matrix_adaptor->store($binding_matrix);

            open my $fh, '<', $motif_dir . '/' . $motif_file;
            my $frequencies_matrix;
            my $line_cnt = 0;

            while ( readline $fh ) {
                my @nucleotide_order = ( 'A', 'C', 'G', 'T' );

                if (/[\d\s]+/) {
                    my @frs = split /\s+/;
                    $frequencies_matrix->{ $nucleotide_order[$line_cnt] }
                        = \@frs;
                    $line_cnt++;
                }
            }

            close $fh;

            my $binding_matrix_frequencies = create_frequencies_objects(
                $logger,     $frequencies_matrix,
                $db_adaptor, $binding_matrix
            );

            my $binding_matrix_frequencies_adaptor = $db_adaptor->get_adaptor(
                'binding_matrix_frequencies_adaptor');

            for my $bmf ( @{binding_matrix_frequencies} ) {
                $binding_matrix_frequencies_adaptor->store($bmf);
            }

        }

    }

}

sub create_matrix_object {
    my ( $logger, $motif_name, $tf_name, $source, $db_adaptor ) = @_;

    my $analysis_adaptor = $db_adaptor->get_adaptor('analysis');
    my $analysis         = $analysis_adaptor->fetch_by_logic_name('SELEX');

    my $feature_type_adaptor = $db_adaptor->get_adaptor('featuretype');
    my $feature_type         = $feature_type_adaptor->fetch_by_name($tf_name);

    if ( !$feature_type ) {
        $logger->warning(
            'Skipping '
                . $motif_name
                . '. Transcription Factor '
                . $tf_name
                . ' not found in database.',
            0, 1
        );
        return;
    }

    my $binding_matrix = Bio::EnsEMBL::Funcgen::BindingMatrix->new(
        'name'         => $motif_name,
        'analysis'     => $analysis,
        'source'       => $source,
        'feature_type' => $feature_type,
    );

    return $binding_matrix;
}

sub create_frequencies_objects {
    my ( $logger, $frequencies_matrix, $db_adaptor, $binding_matrix ) = @_;

    my @binding_matrix_frequencies;
    my $binding_matrix_adaptor = $db_adaptor->get_adaptor('binding_matrix');

    for my $nucleotide ( keys %{$frequencies_matrix} ) {
        for (
            my $position = 0;
            $position < scalar @{ $frequencies_matrix->{$nucleotide} };
            $position++
            )
        {
            my $bmf = Bio::EnsEMBL::Funcgen::BindingMatrixFrequencies->new(
                'binding_matrix' => $binding_matrix,
                'position'       => $position + 1,
                'nucleotide'     => $nucleotide,
                'frequency' =>
                    $frequencies_matrix->{$nucleotide}->[$position],
            );

            push @binding_matrix_frequencies, $bmf;
        }
    }

    return \@binding_matrix_frequencies;
}
