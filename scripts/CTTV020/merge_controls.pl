#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use autodie;
use feature qw(say);

use Data::Printer;

# -----------------------------
# get control (aka input) files
# -----------------------------
my ( %hist_input, %ctcf_input );
open my $fh, '<', $ENV{'STUDY_DIR'} . '/summary.csv';

while ( readline $fh ) {
    if (/hist_input/i) {
        my @fields = split /\t/;
        $fields[0] =~ s/\.fastq\.gz//;    # remove extension
        $hist_input{ $fields[0] } = $fields[10];
    }

    if (/ctcf_input/i) {
        my @fields = split /\t/;
        $fields[0] =~ s/\.fastq\.gz//;    # remove extension
        $ctcf_input{ $fields[0] } = $fields[10];
    }
}

close $fh;

# -----------------------------
# find the ss ones (subsampled)
# -----------------------------
my ( $hist_to_merge, $ctcf_to_merge );

for my $file ( keys %hist_input ) {
    if ( -e $file . '_ss.fastq.gz' ) {
        $hist_to_merge .= $file . '_ss.fastq.gz ';
    }
    elsif ( -e $file . '.fastq.gz' ) {
        $hist_to_merge .= $file . '.fastq.gz ';
    }
    else {
        say STDERR 'File not found: ' . $file;
    }
}

for my $file ( keys %ctcf_input ) {
    if ( -e $file . '_ss.fastq.gz' ) {
        $ctcf_to_merge .= $file . '_ss.fastq.gz ';
    }
    elsif ( -e $file . '.fastq.gz' ) {
        $ctcf_to_merge .= $file . '.fastq.gz ';
    }
    else {
        say STDERR 'File not found: ' . $file;
    }
}

my $hist_merge_cmd
    = 'bsub -e hist_merge.err -o hist_merge.out \'cat ' . $hist_to_merge . ' > hist_input.fastq.gz\'';

my $ctcf_merge_cmd
    = 'bsub -e ctcf_merge.err -o ctcf_merge.out \'cat ' . $ctcf_to_merge . ' > ctcf_input.fastq.gz\'';

my $exit_code = system($hist_merge_cmd);
if ( $exit_code != 0 ) {
    say STDERR 'This system call has failed: ' . $hist_merge_cmd;
}

$exit_code = system($ctcf_merge_cmd);
if ( $exit_code != 0 ) {
    say STDERR 'This system call has failed: ' . $ctcf_merge_cmd;
}
