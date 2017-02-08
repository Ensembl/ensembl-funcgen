#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use autodie;
use feature qw(say);

use Storable;

my $workdir = $ENV{'workdir'};

my @filenames;
my $file_list = shift or die "Please provide a file list";

open my $fh, '<', $file_list;

while (<$fh>) {
    if (/^dataObj: (\S+)\n/) {
        push @filenames, $1;
    }
}

close $fh;

my $data_dir = $ENV{'DATA_DIR'};
store \@filenames, $data_dir . '/filenames_variable';

my $bsub_cmd
    = 'bsub -J download[1-'
    . scalar @filenames
    . ']%15 -o download.%I.out -e download.%I.err \''
    . $workdir
    . '/lib/ensembl-funcgen/scripts/CTTV020/download_job.pl '
    . '\'';
# say $bsub_cmd;
my $exit_code = system($bsub_cmd);
if ( $exit_code != 0 ) {
    say STDERR 'This system call has failed: ' . $bsub_cmd;
}

say 'An array of '
    . scalar @filenames
    . ' download jobs has been submitted to the farm.';

