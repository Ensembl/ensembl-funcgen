#!/usr/bin/env perl

use strict;
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use Carp;
use strict;
use Getopt::Long;

my $species;
my $registry;
my $emissions_file;
my $segmentation;

GetOptions (
   'species=s'          => \$species,
   'registry=s'         => \$registry,
   'emissions_file=s'   => \$emissions_file,
   'segmentation=s'     => \$segmentation,
);

use Bio::EnsEMBL::Utils::Logger;

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

$logger->info("registry          = " . $registry         . "\n");
$logger->info("species           = " . $species          . "\n");
$logger->info("emissions_file    = " . $emissions_file   . "\n");
$logger->info("segmentation      = " . $segmentation     . "\n");

use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->load_all($registry);

my $funcgen_dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
my $segmentation_state_emission_adaptor   = $funcgen_dba->get_SegmentationStateEmissionAdaptor;

open my $fh, '<', $emissions_file or die("Can't open $emissions_file");
my $parsed_emission = parse_emission_matrix($fh);
$fh->close;

foreach my $current_emission (sort { $a->{state} <=> $b->{state} } @$parsed_emission) {

    my $segmentation_state_emisison = Bio::EnsEMBL::Funcgen::SegmentationStateEmission->new(
        -state        => $current_emission->{'state'},
        -segmentation => $segmentation,
        -CTCF         => $current_emission->{'CTCF'},
        -DNase1       => $current_emission->{'DNase1'},
        -H3K27ac      => $current_emission->{'H3K27ac'},
        -H3K27me3     => $current_emission->{'H3K27me3'},
        -H3K36me3     => $current_emission->{'H3K36me3'},
        -H3K4me1      => $current_emission->{'H3K4me1'},
        -H3K4me2      => $current_emission->{'H3K4me2'},
        -H3K4me3      => $current_emission->{'H3K4me3'},
        -H3K9ac       => $current_emission->{'H3K9ac'},
        -H3K9me3      => $current_emission->{'H3K9me3'},
    );
    $segmentation_state_emission_adaptor->store($segmentation_state_emisison);
}

sub parse_emission_matrix {

    my $fh = shift;
    
    my $first_line = <$fh>;
    chomp $first_line;
    my @headers = split "\t", $first_line;
    
    # Change from 'state (Emission order)' in the file.
    $headers[0] = 'state';
    
    my $emission_matrix = parse_tab_separated(
        $fh,
        \@headers,
    );
    
    return $emission_matrix;
}

sub parse_tab_separated {

    my $fh      = shift;
    my $headers = shift;
    my @all;
    
    while (my $current_line = <$fh>) {
    
        chomp $current_line;
        my @values = split "\t", $current_line;
        
        use List::MoreUtils qw( zip );
        my %row = zip(@$headers, @values);
        push @all, \%row;
    }
    return \@all;
}
