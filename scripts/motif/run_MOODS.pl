#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use autodie;
use feature qw(say);

use Data::Printer;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::BindingMatrix;

main();

sub main {
    my $job_index = $ENV{'LSB_JOBINDEX'};

    my $db_adaptor = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(

        -host   => 'mysql-ens-reg-prod-1.ebi.ac.uk',
        -user   => 'ensadmin',
        -pass   => 'ensembl',
        -port   => 4526,
        -dbname => 'ilavidas_motif_homo_sapiens_funcgen_89_38'
    );

    my $binding_matrix_adaptor = $db_adaptor->get_adaptor('BindingMatrix');

    # my $binding_matrices = $binding_matrix_adaptor->fetch_all();
    my $binding_matrix = $binding_matrix_adaptor->fetch_by_dbID($job_index);

    run_MOODS($binding_matrix);

}

sub run_MOODS {
    my ($binding_matrix) = @_;

    my $workdir = '/hps/nobackup/production/ensembl/ilavidas/motif/';

    my $motif_path
        = $workdir . 'matrices/' . $binding_matrix->name() . '.pfm';
    my $target_path = $workdir . 'GRCh38.fa';
    my $output_path
        = $workdir . 'moods_output/' . $binding_matrix->name() . '.out';

    open my $out, '>', $motif_path;

    print $out $binding_matrix->get_frequencies_as_string();

    close $out;

# moods_dna.py -m <motif_name> -s <target_fasta> -p 0.01 >output_file
    my $command
        = 'moods_dna.py -m '
        . $motif_path . ' -s '
        . $target_path . ' -p '
        . $binding_matrix->threshold() . ' >'
        . $output_path;

    print $command . "\n";

    my $exit_code = system($command);

    if ($exit_code) {
        say 'COULD NOT EXECUTE ' . $command;
    }
}
