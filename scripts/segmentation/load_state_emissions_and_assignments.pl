#!/usr/bin/env perl

use strict;
use Bio::EnsEMBL::Registry;
use v5.10;
use Data::Dumper;
use Carp;
$Data::Dumper::Maxdepth = 3;

=head1 

/homes/mnuhn/various_scripts/new_ersa/load_state_emissions.pl


drop table if exists segmentation_state_emission;

create table segmentation_state_emission (
    segmentation_state_emission_id int(27) unsigned NOT NULL AUTO_INCREMENT,
    state int(7) NOT NULL,
    CTCF double NOT NULL,
    DNase1 double NOT NULL,
    H3K27ac double NOT NULL,
    H3K27me3 double NOT NULL,
    H3K36me3 double NOT NULL,
    H3K4me1 double NOT NULL,
    H3K4me2 double NOT NULL,
    H3K4me3 double NOT NULL,
    H3K9ac double NOT NULL,
    H3K9me3 double NOT NULL,
    PRIMARY KEY (segmentation_state_emission_id)
);

drop table if exists segmentation_state_assignment;

create table segmentation_state_assignment (
    segmentation_state_assignment_id int(35) unsigned NOT NULL AUTO_INCREMENT,
    state int(7) NOT NULL,
    segmentation varchar(100) NOT NULL,
    assignment varchar(100) NOT NULL,
    PRIMARY KEY (segmentation_state_assignment_id)
);

emissions_file=/gpfs/nobackup/ensembl/mnuhn/mnuhn/regulatory_build_pipeline_run3/temp_dir/segmentation/mus_musculus/learn_model/emissions_25.txt
assignments_file=/gpfs/nobackup/ensembl/mnuhn/mnuhn/regulatory_build_pipeline_run3/temp_dir/regulatory_build/mus_musculus/tmp/assignments.txt

perl scripts/segmentation//load_state_emissions_and_assignments.pl \
    --species          mus_musculus \
    --registry         /homes/mnuhn/work_dir_ersa/lib/ensembl-funcgen/registry.pm \
    --emissions_file   $emissions_file \
    --assignments_file $assignments_file

perl scripts/segmentation/segmentation_statistics.pl \
    --species          mus_musculus \
    --registry         /homes/mnuhn/work_dir_ersa/lib/ensembl-funcgen/registry.pm \
    --output_directory ./segmentation_report

=cut 

use strict;
use Getopt::Long;

my $species;
my $registry;
my $emissions_file;
my $segmentation;
my $assignments_file;

GetOptions (
   'species=s'          => \$species,
   'registry=s'         => \$registry,
   'emissions_file=s'   => \$emissions_file,
   'segmentation=s'     => \$segmentation,
   'assignments_file=s' => \$assignments_file,
);

use Bio::EnsEMBL::Utils::Logger;

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

$logger->info("registry          = " . $registry         . "\n");
$logger->info("species           = " . $species          . "\n");
$logger->info("emissions_file    = " . $emissions_file   . "\n");
$logger->info("segmentation      = " . $segmentation     . "\n");
$logger->info("assignments_file  = " . $assignments_file . "\n");

use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->load_all($registry);

my $funcgen_dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');

my $dbc = $funcgen_dba->dbc;

# $dbc->do("truncate segmentation_state_assignment;");
# $dbc->do("truncate segmentation_state_emission;");

my $segmentation_state_emission_adaptor   = $funcgen_dba->get_SegmentationStateEmissionAdaptor;
my $segmentation_state_assignment_adaptor = $funcgen_dba->get_SegmentationStateAssignmentAdaptor;

open my $fh, '<', $emissions_file or die("Can't open $emissions_file");
my $parsed_emission = parse_emission_matrix($fh);
$fh->close;

# print Dumper($parsed);

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

if (! defined $assignments_file) {
    exit(0);
}

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
