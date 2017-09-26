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
use Bio::EnsEMBL::Funcgen::DBSQL::TranscriptionFactorComplexAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::BindingMatrixAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::BindingMatrixFrequenciesAdaptor;

use Bio::EnsEMBL::Funcgen::TranscriptionFactor;
use Bio::EnsEMBL::Funcgen::TranscriptionFactorComplex;
use Bio::EnsEMBL::Funcgen::BindingMatrix;
use Bio::EnsEMBL::Funcgen::BindingMatrixFrequencies;

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

        my $tf1_name     = $fields[0];
        my $tf2_name     = $fields[1];
        my $gene_id1     = $fields[2];
        my $gene_id2     = $fields[3];
        my @matrix_names = split /,/, $fields[8];    #representative motifs

        # if($tf1_name eq 'ARX'){die;}

        # -------------------------
        # store to funcgen database
        # -------------------------
        my $first_transcription_factor
            = store_transcription_factor( $tf1_name, $gene_id1, $db_adaptor,
            $logger );
        my $second_transcription_factor
            = store_transcription_factor( $tf2_name, $gene_id2, $db_adaptor,
            $logger );

        my $transcription_factor_complex = store_transcription_factor_complex( $first_transcription_factor,
            $second_transcription_factor, $db_adaptor, $logger );

        for my $matrix_name (@matrix_names){
            store_binding_matrix( $matrix_name, $source, $db_adaptor, $logger, $motif_dir, $transcription_factor_complex);
        }
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

    if ( !$transcription_factor_complex ) {

        $transcription_factor_complex
            = Bio::EnsEMBL::Funcgen::TranscriptionFactorComplex->new(
            -PRODUCTION_NAME => $production_name,
            -DISPLAY_NAME    => $display_name,
            -COMPONENTS      => \@components,
            );

        $transcription_factor_complex_adaptor->store(
            $transcription_factor_complex);

        $logger->info(
            'Transcription Factor Complex '
                . $display_name
                . ' stored' . "\n",
            0, 1
        );

    }
    else {
        $logger->info(
            'Transcription Factor Complex '
                . $display_name
                . ' found in DB. Skipping...' . "\n",
            0, 1
        );
    }

    return $transcription_factor_complex;
}


sub store_binding_matrix {
    my ( $matrix_name, $source, $db_adaptor, $logger, $motif_dir,
        $transcription_factor_complex )
        = @_;

    my $binding_matrix_adaptor = $db_adaptor->get_adaptor('BindingMatrix');

    my $binding_matrix = $binding_matrix_adaptor->fetch_by_name($matrix_name);

    if ( !$binding_matrix ) {

        my $frequencies
            = read_matrix_file( $matrix_name, $logger, $motif_dir );
        if ( !$frequencies ) { return; }

        $binding_matrix = Bio::EnsEMBL::Funcgen::BindingMatrix->new(
            -NAME        => $matrix_name,
            -THRESHOLD   => 0.1,
            -SOURCE      => $source,
            -FREQUENCIES => $frequencies,
            -ASSOCIATED_TRANSCRIPTION_FACTORS =>
                [$transcription_factor_complex],
        );

        $binding_matrix_adaptor->store($binding_matrix);

        $logger->info( 'Binding Matrix ' . $matrix_name . ' stored' . "\n",
            0, 1 );
    }
    else {
        $binding_matrix->{associated_transcription_factors} = [$transcription_factor_complex ];
        $binding_matrix_adaptor->_store_binding_matrix_transcription_factor_complex($binding_matrix);
        $logger->info(
            'Binding Matrix '
                . $matrix_name
                . ' found in DB. Skipping...' . "\n",
            0, 1
        );
    }
}


sub read_matrix_file {
    my ( $matrix_name, $logger, $motif_dir ) = @_;
    my $line_cnt=0;
    my @frequencies_matrix;

    my $filepath = $motif_dir . '/' . $matrix_name . '.pfm';
    
    if (! -e $filepath){
        $logger->info($filepath . ' not found' . "\n", 0,1);
        return;
    }

    open my $fh, '<', $filepath;

    while ( readline $fh ) {
        my @nucleotide_order = ( 'A', 'C', 'G', 'T' );

        if (/[\d\s]+/) {
            my @frs = split /\s+/;

            my $col_cnt = 1;

            for my $fr (@frs) {
                $frequencies_matrix[$col_cnt] //= {};
                $frequencies_matrix[$col_cnt]
                    ->{ $nucleotide_order[$line_cnt] } = $fr;
                $col_cnt++;
            }

            $line_cnt++;
        }
    }

    close $fh;

    return \@frequencies_matrix;

}