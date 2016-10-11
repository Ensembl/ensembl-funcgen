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
        exit 1;
    }

    my $json_dir = abs_path( $ARGV[0] );
    my $data     = get_data_from_jsons($json_dir);

    while ( my ( $json_file, $json_data ) = each %{$data} ) {
        check_mandatory_fields( $json_file, $json_data );
    }

    print_csv($data);

    say 'success';
}

sub print_csv {
    my ( $data, $output ) = @_;

    for my $j ( values %{$data} ) {
        print $j->{accession} . ","
            . $j->{replicate}->{experiment}->{biosample_term_name} . ","
            . $j->{replicate}->{experiment}->{target}->{label} . ","
            . $j->{replicate}->{biological_replicate_number} . ","
            . $j->{replicate}->{technical_replicate_number} . ","
            . $j->{md5sum} . ","
            . $j->{funcgen}->{local_url} . ","
            . $j->{funcgen}->{is_control} . ","
            . $j->{replicate}->{experiment}->{assay_title} . ","
            . $j->{funcgen}->{experimental_group} . ",";

        print join( ',',
            @{ $j->{replicate}->{experiment}->{biosample_term_id} } );
        print ",";

        if ( $j->{funcgen}->{is_control} == 0 ) {
            print join( ',', @{ $j->{controlled_by} } );
            print ",";
        }

        print join( ',',
            @{ $j->{replicate}->{experiment}->{biosample_synonyms} } );

        print "\n";
    }

}

sub get_data_from_jsons {
    my ($json_dir) = @_;

    my %data;

    opendir( my $dir_h, $json_dir );
    my @json_files = grep {/.json/} readdir $dir_h;

    if ( scalar @json_files == 0 ) {
        say STDERR "ERROR: No json files found in $json_dir";
        exit 1;
    }

    for my $json_file (@json_files) {
        $data{$json_file} = json_file_to_perl( $json_dir . '/' . $json_file );
    }

    return \%data;
}

sub check_mandatory_fields {
    my ( $json_file, $json_data ) = @_;
    my $multi_error_flag = 0;

    my @mandatory_fields = (
        "accession",
        "md5sum",
        "funcgen.experimental_group",
        "funcgen.local_url",
        "funcgen.is_control",
        "replicate.biological_replicate_number",
        "replicate.experiment.assay_title",
        "replicate.experiment.biosample_term_name",
        "replicate.experiment.biosample_term_id",
        "replicate.experiment.target.label",
        "replicate.technical_replicate_number",
    );

    for my $man_field (@mandatory_fields) {
        my $error_flag = 0;
        my @tiers = split /\./, $man_field;

        my $field = $json_data;
        for ( my $i = 0; $i <= $#tiers; $i++ ) {
            $field = $field->{ $tiers[$i] };

            if ( !defined $field || length($field) == 0 ) {
                $error_flag       = 1;
                $multi_error_flag = 1;
            }
        }

        say STDERR
            "ERROR: In file '$json_file' mandatory json attribute '$man_field' is not defined!"
            if $error_flag == 1;
    }

    if ($multi_error_flag) {
        say STDERR "Exiting... Please check this file: $json_file";
        exit 1;
    }

}

sub usage {
    say STDERR "Usage: perl json2csv.pl <json_dir>";
    say STDERR "Please set the directory which contains the json files";
}
