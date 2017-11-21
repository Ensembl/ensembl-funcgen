#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use autodie;
use feature qw(say);

use Data::Printer;
use File::Basename;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Funcgen::DBSQL::MotifFeatureAdaptor;
use Bio::EnsEMBL::Funcgen::MotifFeature;

main();

sub main {
    my $job_index = $ENV{'LSB_JOBINDEX'};

    my $funcgen_adaptor = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
        -host   => 'mysql-ens-reg-prod-1.ebi.ac.uk',
        -user   => 'ensadmin',
        -pass   => 'ensembl',
        -port   => 4526,
        -dbname => 'ilavidas_motif_2_homo_sapiens_funcgen_89_38'
    );

    my $core_adaptor = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -species => 'homo_sapiens',
        -group   => 'core',
        -host    => 'ensembldb.ensembl.org',
        -user    => 'anonymous',
        -dbname  => 'homo_sapiens_core_90_38',
    );

    my $moods_output_dir = '/hps/nobackup/production/ensembl/ilavidas/motif/moods_output/';
    opendir (my $dirh, $moods_output_dir);
    my @files = grep {!/^\./} readdir $dirh;

    my $moods_output_file = $files[$job_index-1];


    my ($binding_matrix_name) = fileparse($moods_output_file);
    $binding_matrix_name =~ s/\.out//g;

    my $binding_matrix_adaptor
        = $funcgen_adaptor->get_adaptor('BindingMatrix');

    my $binding_matrix
        = $binding_matrix_adaptor->fetch_by_name($binding_matrix_name);

    my $matrix_length = $binding_matrix->length();

    my $threshold = $binding_matrix->threshold();

    my $slice_adaptor = $core_adaptor->get_adaptor('Slice');

    my $motif_feature_adaptor = $funcgen_adaptor->get_adaptor('MotifFeature');

    open my $fh, '<', $moods_output_file;

    while ( readline $fh ) {
        chomp;

        my @fields = split /,/;

        my @fasta_header = split /\s+/, $fields[0];

        my $score = $fields[4];

        # skip if score is below the defined matrix threshold
        if ( $score < $threshold ) { next; }

        my $seq_region_name = pop @fasta_header;

        # MOODS first position of sequence is position 0
        my $seq_region_start = $fields[2] + 1;
        my $seq_region_end   = $seq_region_start + $matrix_length - 1;

        my $seq_region_strand = $fields[3];
        if ( $seq_region_strand eq '+' ) {
            $seq_region_strand = 1;
        }
        elsif ( $seq_region_strand eq '-' ) {
            $seq_region_strand = -1;
        }
        else {
            die 'Invalid seq_region_strand: ' . $seq_region_strand;
        }

        my $moods_sequence = $fields[5];

        my $slice = $slice_adaptor->fetch_by_name($seq_region_name);

        my $motif_feature = Bio::EnsEMBL::Funcgen::MotifFeature->new(
            -SLICE          => $slice,
            -START          => $seq_region_start,
            -END            => $seq_region_end,
            -STRAND         => $seq_region_strand,
            -BINDING_MATRIX => $binding_matrix,
            -SCORE          => $score,
            -STABLE_ID      => 'ENSM0000001'
        );

        # quick on-the-fly healthcheck
        # make sure that the slice sequence is identical
        # to the one reported by MOODS
        my $slice_sequence
            = $slice->subseq( $seq_region_start, $seq_region_end );

        if ( $slice_sequence ne $moods_sequence ) {
            say STDERR 'Sequences do not match:';
            say STDERR 'MOODS: ' . $moods_sequence;
            say STDERR 'SLICE: ' . $slice_sequence;
        }

        $motif_feature_adaptor->store($motif_feature);
    }
}
