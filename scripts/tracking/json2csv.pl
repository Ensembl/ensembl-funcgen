#!/usr/bin/env perl

use strict;
use warnings;

# use diagnostics;
use autodie;
use feature qw(say);

use Data::Printer;
use Cwd 'abs_path';

# use JSON qw(decode_json);
# use Path::Tiny;
use JSON::Parse 'json_file_to_perl';

main();

sub main {

    if ( !$ARGV[0] ) {
        usage();
        exit 0;
    }

    my $json_dir = abs_path( $ARGV[0] );
    my $data     = get_data_from_jsons($json_dir);

    p $data->[0]->{replicate}->{experiment}->{files}->[1];

}

sub get_data_from_jsons {
    my ($json_dir) = @_;
    my @data;

    opendir( my $dir_h, $json_dir );
    my @json_files = grep {/.json/} readdir $dir_h;

    if ( scalar @json_files == 0 ) {
        say "No json files found in $json_dir";
        exit 0;
    }

    for my $json_file (@json_files) {
        push @data, json_file_to_perl( $json_dir . '/' . $json_file );
    }

    return \@data;
}

sub check_mandatory_fields {

}

sub usage {
    say "Usage: perl json2csv.pl <json_dir>";
    say "Please set the directory which contains the json files";
}
