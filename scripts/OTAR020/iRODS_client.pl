#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use autodie;
use feature qw(say);

# use Data::Printer;
# use Storable;
use Expect;
use Term::ReadKey;

#TODO set study as a parameter, no need to check all studies every time

my $workdir = '/lustre/scratch117/ensembl/il4/OTAR020/';

$ENV{'SANG_PASS'} = get_password();

my %OTAR020_dirs = (

    'CTTV020 epigenomes of cell lines PILOT'       => $workdir . 'PILOT',
    'CTTV020 epigenomes of cell lines PILOT - RNA' => $workdir . 'PILOT_RNA',

    'CTTV020 epigenomes of cell lines FULL - ChIP' => $workdir . 'FULL_ChIP',
    'CTTV020 epigenomes of cell lines FULL - RNA'  => $workdir . 'FULL_RNA',
    'CTTV020 Epigenomes of Cell Lines ATAC-seq'    => $workdir . 'ATAC',

    'CTTV020 Epigenomes of Cell Lines ChIPmentation' => $workdir
        . 'ChIPmentation',
);

main();

sub main {

    # kinit();

    for my $study ( keys %OTAR020_dirs ) {

        my @new_datasets;
        my @crams;

        # query iRODS server for all cram files related to this study
        my $stdout
            = `imeta qu -z seq -d study = "$study" and target = 1 and manual_qc = 1`;

        # read stdout line by line
        foreach ( split( /\n/, $stdout ) ) {

            # identify cram filenames and store them in an array
            if (/\s+(\S+\.cram)/) {
                push @crams, $1;
            }
        }

        for my $cram (@crams) {

           # if the cram file has not already been downloaded mark it as "new"
            unless ( -e $OTAR020_dirs{$study} . '/cram/' . $cram ) {
                push @new_datasets, $cram;
            }
        }

        if ( scalar @new_datasets > 0 ) {

            print 'There are '
                . scalar @new_datasets
                . ' new datasets for "'
                . $study
                . '". Generate meta files (y/N)? ';

            my $stdin = readline();
            chomp $stdin;

            if ( $stdin eq 'y' ) {
                generate_meta( $study, \@new_datasets );
            }

            print 'There are '
                . scalar @new_datasets
                . ' new datasets for "'
                . $study
                . '". Download (y/N)? ';

            $stdin = readline();
            chomp $stdin;

            if ( $stdin eq 'y' ) {
                submit_download_array_job( $study, \@new_datasets );
            }

        }
        else {
            say 'There are no new datasets for "' . $study . '"';
        }
    }
}

sub generate_meta {
    my ( $study, $datasets ) = @_;

    for my $dataset ( @{$datasets} ) {
        my ($parent_dir) = split /_/, $dataset;

        my $command
            = 'imeta ls -d /seq/'
            . $parent_dir . '/'
            . $dataset . ' >'
            . $OTAR020_dirs{$study}
            . '/meta/'
            . $dataset . '.meta';

        say 'Generating meta file for ' . $dataset;

        my $exit_code = system($command);
        if ( $exit_code != 0 ) {
            say 'This system call has failed: ' . $command;
        }
    }
}

sub submit_download_array_job {
    my ( $study, $datasets ) = @_;

    my $datasets_string = join( ",", @{$datasets} );

    my $bsub_cmd
        = 'bsub -J download[1-'
        . scalar @{$datasets}
        . ']%15 -o '
        . $OTAR020_dirs{$study}
        . '/download.%I.out -e '
        . $OTAR020_dirs{$study}
        . '/download.%I.err \'' . 'perl '
        . $workdir
        . '/lib/ensembl-funcgen/scripts/OTAR020/download_job.pl ' . '"'
        . $datasets_string . '" '
        . $OTAR020_dirs{$study}
        . '/cram/' . '\'';

    my $exit_code = system($bsub_cmd);
    if ( $exit_code != 0 ) {
        say STDERR 'This system call has failed: ' . $bsub_cmd;
    }

    say 'An array of '
        . scalar @{$datasets}
        . ' download jobs has been submitted to the farm.';
}

sub get_password {
    print 'Please type your Sanger password:';

    ReadMode('noecho');    # make password invisible on terminal
    my $password = ReadLine(0);
    chomp $password;
    ReadMode(0);           # restore typing visibility on terminal
    print "\n";

    return $password;
}

sub kinit {
    my $session = Expect->spawn('kinit') or die $!;
    if ($session->expect(
            5, 'Password for ' . $ENV{'USER'} . '\@INTERNAL.SANGER.AC.UK:'
        )
        )
    {

        print $session $ENV{'SANG_PASS'} . "\r";
    }
    else {
        say STDERR 'kinit failed';
    }

    $session->soft_close();
}
