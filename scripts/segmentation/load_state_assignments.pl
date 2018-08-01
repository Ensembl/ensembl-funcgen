#!/usr/bin/env perl

use strict;
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use Carp;
use strict;
use Getopt::Long;

my $species;
my $registry;
my $assignments_file;

GetOptions (
   'species=s'          => \$species,
   'registry=s'         => \$registry,
   'assignments_file=s' => \$assignments_file,
);

use Bio::EnsEMBL::Utils::Logger;

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

$logger->info("registry          = " . $registry         . "\n");
$logger->info("species           = " . $species          . "\n");
$logger->info("assignments_file  = " . $assignments_file . "\n");

use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->load_all($registry);

my $funcgen_dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
my $segmentation_state_assignment_adaptor = $funcgen_dba->get_SegmentationStateAssignmentAdaptor;

open $fh, '<', $assignments_file || die("Can't open $assignments_file");

my $parsed_assignment = parse_assignments(
    $fh,
    [ 'segmentation', 'e_state', 'assignment' ]
);
$fh->close;

foreach my $current_assignment (sort { $a->{state} <=> $b->{state} } @$parsed_assignment) {
    
    my $segmentation_state_assignment = Bio::EnsEMBL::Funcgen::SegmentationStateAssignment->new(
        -state        => $current_assignment->{'state'},
        -segmentation => $current_assignment->{'segmentation'},
        -assignment   => $current_assignment->{'assignment'},
    );
    $segmentation_state_assignment_adaptor->store($segmentation_state_assignment);

}

sub parse_assignments {

    my $fh = shift;
    my $headers = shift;
    
    my $assignments = parse_tab_separated(
        $fh,
        $headers,
    );
    map {
        my $row = $_;
        my $e_state = $row->{e_state};
        my $state_found = $e_state =~ /^E(.+$)/;
        if (! $state_found) {
            die("Couldn't find state number in $_!\n" . Dumper($row));
        }
        my $state = $1;
        $row->{'state'} = $state;
    } @$assignments;
    
    return $assignments;
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
