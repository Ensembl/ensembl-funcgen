#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use autodie;
use feature qw(say);
use Data::Printer;

# ----------------------------------------------------------------
# Open every meta file and get the attributes of interest.
# ----------------------------------------------------------------
opendir my $meta_dirh, $ARGV[0];
my @meta_files = grep { !/^\./ } readdir $meta_dirh;    # exclude hidden files

my %data;
my @attributes_of_interest = (
    'sample_public_name', 'id_run',
    'lane',               'tag_index',
    'sample_id',          'total_reads'
);

foreach my $meta_file (@meta_files) {

    open my $mfh, '<', $ARGV[0] . '/' . $meta_file;
    $/ = "----\n";    # read entry by entry, not line by line

    while (<$mfh>) {
        foreach my $attribute (@attributes_of_interest) {
            if (/attribute: $attribute\nvalue: (.*)\nunits/) {
                $data{$meta_file}->{$attribute} = $1;
            }
        }
    }

    close $mfh;

}

# -----------------
# Get md5 checksums
# -----------------
open my $md5_fh, '<', $ENV{'WAREHOUSE_DIR'} . '/md5.txt';
$/ = "\n";    # read line by line

my %md5;

while ( readline $md5_fh ) {
    if (/^(\S+)\s+(\S+)/) {
        $md5{$2} = $1;
    }
}

close $md5_fh;

# --------------
# Print csv file
# --------------
open my $out, '>', $ENV{'STUDY_DIR'} . '/summary.csv';

foreach my $key ( sort keys %data ) {
    ( my $file = $key ) =~ s/meta/fastq\.gz/;

    print $out $file . "\t";

    foreach my $attribute (@attributes_of_interest) {

        if ( defined $data{$key}->{$attribute} ) {
            if ( $attribute eq 'sample_public_name' ) {
                break_sample_name( $data{$key}->{$attribute} );
            }
            else {
                print $out $data{$key}->{$attribute} . "\t";
            }
        }
        else {
            print $out 'NULL' . "\t";
        }
    }
    print $out $md5{$file} . "\t";
    print $out "\n";
}

close $out;

sub break_sample_name {
    my ($sample_public_name) = @_;

    my ( $cell_line, $br, $tr, $feature ) = split /_/, $sample_public_name;

    print $out $cell_line . "\t" . $br . "\t" . $tr . "\t" . $feature . "\t";
}

