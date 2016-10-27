#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use autodie;
use feature qw(say);

use Storable;

my $study_dir = $ENV{'STUDY_DIR'};
my $job_index = $ENV{'LSB_JOBINDEX'};

my @filenames = @{ retrieve( $study_dir . '/filenames_variable' ) };

my ( $filename, $extension ) = split /\./, $filenames[ $job_index - 1 ];

# my $command
#     = 'cat '
#     . $filenames[ $job_index - 1 ]
#     . ' | /software/solexa/pkg/biobambam/current/bin/bamtofastq inputformat=cram exclude=SECONDARY,SUPPLEMENTARY,QCFAIL'
#     . ' | gzip >'
#     . $filename
#     . '.fastq.gz';
my $command
    = '/software/ensembl/funcgen//samtools bam2fq '
    . $filenames[ $job_index - 1 ]
    . ' | gzip >'
    . $filename
    . '.fastq.gz';

my $exit_code = system($command);
if ( $exit_code != 0 ) {
    say STDERR 'This system call has failed: ' . $command;
}

my $md5_cmd = 'md5sum ' . $filename . '.fastq.gz >>md5.txt';

$exit_code = system($md5_cmd);
if ( $exit_code != 0 ) {
    say STDERR 'This system call has failed: ' . $md5_cmd;
}
