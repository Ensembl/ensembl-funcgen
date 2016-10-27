#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use autodie;
use feature qw(say);

use Data::Printer;
use Storable;

my $cttv020_dir = $ENV{'CTTV020_DIR'};
my $study_dir   = $ENV{'STUDY_DIR'};

my @filenames = @{ retrieve( $study_dir . '/filenames_variable' ) };
# s/\.cram// for @filenames;    # remove .cram extension from filenames

my $bsub_cmd
    = 'bsub -M6000 -R"select[mem>6000] rusage[mem=6000]" -J convert[1-'
    . scalar @filenames
    . '] -o convert.%I.out -e convert.%I.err \''
    . $cttv020_dir
    . '/lib/ensembl-funcgen/scripts/CTTV020/convertCram2fq_job.pl ' . '\'';

my $exit_code = system($bsub_cmd);
if ( $exit_code != 0 ) {
    say STDERR 'This system call has failed: ' . $bsub_cmd;
}

say 'An array of '
    . scalar @filenames
    . ' conversion jobs has been submitted to the farm.';

