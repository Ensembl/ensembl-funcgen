#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use autodie;
use feature qw(say);

use Storable;

my $cttv020_dir = $ENV{'CTTV020_DIR'};

my @filenames;
my $file_list = shift or die "Please provide a file list";

open my $fh, '<', $file_list;

while (<$fh>) {
    if (/^dataObj: (\S+)\n/) {
        push @filenames, $1;
    }
}

close $fh;

my $study_dir = $ENV{'STUDY_DIR'};
store \@filenames, $study_dir . '/filenames_variable';

my $bsub_cmd
    = 'bsub -J download[1-'
    . scalar @filenames
    . ']%15 -o download.%I.out -e download.%I.err \''
    . $cttv020_dir
    . '/scripts/download_job.pl '
    . '\'';
# say $bsub_cmd;
my $exit_code = system($bsub_cmd);
if ( $exit_code != 0 ) {
    say STDERR 'This system call has failed: ' . $bsub_cmd;
}

say 'An array of '
    . scalar @filenames
    . ' download jobs has been submitted to the farm.';

