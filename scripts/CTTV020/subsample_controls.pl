#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use autodie;
use feature qw(say);
use Data::Printer;

use List::Util qw( min );
use Math::Round qw/round/;

# -----------------------------
# get control (aka input) files
# -----------------------------
my ( %his_input, %ctcf_input );
open my $fh, '<', $ENV{'STUDY_DIR'} . '/summary.csv';

while ( readline $fh ) {
    if (/his_input/i) {
        my @fields = split /\t/;
        $fields[0] =~ s/\.fastq\.gz//;    # remove extension
        $his_input{ $fields[0] } = $fields[10];
    }

    if (/ctcf_input/i) {
        my @fields = split /\t/;
        $fields[0] =~ s/\.fastq\.gz//;    # remove extension
        $ctcf_input{ $fields[0] } = $fields[10];
    }
}

close $fh;

# -----------------------------
# calculate ratio for each file
# -----------------------------
my ( %his_ratio, %ctcf_ratio );

my $his_min  = min values %his_input;
my $ctcf_min = min values %ctcf_input;

for my $key ( keys %his_input ) {
    my $ratio = $his_min / $his_input{$key} * 100;
    $ratio = round($ratio) / 100;
    $his_ratio{$key} = $ratio;
}

for my $key ( keys %ctcf_input ) {
    my $ratio = $ctcf_min / $ctcf_input{$key} * 100;
    $ratio = round($ratio) / 100;
    $ctcf_ratio{$key} = $ratio;
}

# p %his_ratio;
# p %ctcf_ratio;

# p $ctcf_min;

# ------------------------------
# submit subsampling job to farm
# ------------------------------
# samtools view -s 0.96 -b 16723_2#28.cram | samtools bam2fq - | gzip > 16723_2#28_96p.fastq.gz
my $cnt;
for my $file ( keys %his_ratio ) {
    if ( $his_ratio{$file} != 1 ) {

        $cnt++;

        my $subsample_cmd
            = 'samtools view -s '
            . $his_ratio{$file} . ' -b '
            . $file
            . '.cram | samtools bam2fq - | gzip >'
            . $file
            . '_ss.fastq.gz';

        my $bsub_cmd
            = "bsub -M6000 -R\"select[mem>6000] rusage[mem=6000]\" -e subsample.$cnt.err -o subsample.$cnt.out \'$subsample_cmd\'";

        my $exit_code = system($bsub_cmd);
        if ( $exit_code != 0 ) {
            say STDERR 'This system call has failed: ' . $bsub_cmd;
        }
    }
}

for my $file ( keys %ctcf_ratio ) {
    if ( $ctcf_ratio{$file} != 1 ) {

        $cnt++;

        my $subsample_cmd
            = 'samtools view -s '
            . $ctcf_ratio{$file} . ' -b '
            . $file
            . '.cram | samtools bam2fq - | gzip >'
            . $file
            . '_ss.fastq.gz';

        my $bsub_cmd
            = "bsub -M6000 -R\"select[mem>6000] rusage[mem=6000]\" -e subsample.$cnt.err -o subsample.$cnt.out \'$subsample_cmd\'";

        my $exit_code = system($bsub_cmd);
        if ( $exit_code != 0 ) {
            say STDERR 'This system call has failed: ' . $bsub_cmd;
        }
    }
}

say $cnt . ' jobs have been submitted to the farm';
