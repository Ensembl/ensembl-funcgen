#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use autodie;
use feature qw(say);

use Data::Printer;

# -----------------
# Get md5 checksums
# -----------------
open my $md5_fh, '<', $ENV{'WAREHOUSE_DIR'} . '/md5.txt';

my %md5;

while ( readline $md5_fh ) {
    if (/^(\S+)\s+(\S+)/) {
        $md5{$2} = $1;
    }
}

close $md5_fh;

# -------------------------------------
# Print non control registration list
# -------------------------------------
my @ct_names;    # cell type names

open my $nonControlOut, '>',
    $ENV{'STUDY_DIR'} . '/non_control_registration_list.csv';
open my $controlOut, '>',
    $ENV{'STUDY_DIR'} . '/control_registration_list.csv';

open my $fh, '<', $ENV{'STUDY_DIR'} . '/summary.csv';

my $header = readline $fh;    # remove header

while ( readline $fh ) {
    my @fields = split /\t/;

    my $filename     = $fields[0];
    my $cell_line    = $fields[1];
    my $feature_type = $fields[2];
    my $br           = $fields[3];
    my $fg_rep       = $fields[5];    # CAUTION: do not use the tr column

    my $cell_type_name
        = create_cell_type_name( $cell_line, $feature_type, $br );

    push @ct_names, $cell_type_name
        unless grep { $_ eq $cell_type_name } @ct_names;

    if ( $feature_type !~ /input/ ) {
        say $nonControlOut $filename . "\t"
            . $cell_type_name . "\t"
            . $feature_type . "\t"
            . $fg_rep . "\t"
            . $md5{$filename};
    }

}

# -------------------------------------
# Print control registration list
# -------------------------------------
# We have to assign a control file to each cell_type

for my $ct_name (@ct_names) {
    my $filename
        = $ct_name =~ /\:TF\:/
        ? 'ctcf_input.fastq.gz'
        : 'hist_input.fastq.gz';

    say $controlOut $filename . "\t"
        . $ct_name . "\t" . 'WCE' . "\t" . '1' . "\t"
        . $md5{$filename};
}

close $nonControlOut;
close $controlOut;
close $fh;

# -------------------------------------
# Subroutine for compiling the ct name
# -------------------------------------

sub create_cell_type_name {
    my ( $cell_line, $feature_type, $br ) = @_;

    my $cell_type_name;

    if ( $feature_type =~ /CTCF/ ) {
        $cell_type_name = $cell_line . ':TF:' . $br;
    }
    elsif ( $feature_type =~ /^H\d/ || $feature_type =~ /^hist/ ) {
        $cell_type_name = $cell_line . ':hist:' . $br;
    }
    else {
        say STDERR 'Feature type not recognized!';
    }

    return $cell_type_name;
}
