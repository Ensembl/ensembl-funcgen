#!/usr/bin/env perl

use strict;
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use Carp;
use Bio::EnsEMBL::Utils::Logger;

=head1 

time mysql $(r2-w details mysql) mnuhn_testdb5_mus_musculus_funcgen_91_38 -e "drop table peak_calling_statistic"

time mysql $(r2-w details mysql) mnuhn_testdb5_mus_musculus_funcgen_91_38 -e "
    create table 
        peak_calling_statistic as
    select 
        peak_calling_id, 
        sum(peak.seq_region_end - peak.seq_region_start + 1) as total_length,
        count(peak.peak_id) as num_peaks,
        avg(peak.seq_region_end - peak.seq_region_start + 1) as average_length
    from 
        peak_calling 
        left join peak using (peak_calling_id) 
    group by 
        peak_calling_id
    ;
"

peak_calling_statistics.pl \
    --species          mus_musculus \
    --registry         /homes/mnuhn/work_dir_ersa/lib/ensembl-funcgen/registry.pm \
    --output_directory /hps/nobackup/production/ensembl/mnuhn/regulatory_build_pipeline_run2/reports/mus_musculus

=cut

use strict;


use Getopt::Long;

my $species;
my $registry;
my $output_directory;

GetOptions (
   'species=s'          => \$species,
   'registry=s'         => \$registry,
   'output_directory=s' => \$output_directory,
);

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

$logger->info("registry          = " . $registry         . "\n");
$logger->info("species           = " . $species          . "\n");
$logger->info("output_directory  = " . $output_directory . "\n");

use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->load_all($registry);

my $mouse_funcgen_dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');

my $peak_calling_statistic_adaptor = $mouse_funcgen_dba->get_PeakCallingStatisticAdaptor;

my $epigenome_adaptor    = $mouse_funcgen_dba->get_EpigenomeAdaptor;
my $feature_type_adaptor = $mouse_funcgen_dba->get_FeatureTypeAdaptor;
my $peak_calling_adaptor = $mouse_funcgen_dba->get_PeakCallingAdaptor;

my $peak_calling_statistics = $peak_calling_statistic_adaptor->fetch_all;

my $file = __FILE__;
use File::Basename qw( dirname basename );
my $description_template = dirname($file) . '/../../templates/peak_calling/report.html';

if (! -e $description_template) {
    die("Can't find $description_template");
}

use File::Path qw( make_path );
make_path( $output_directory );

use Template;
my $tt = Template->new( ABSOLUTE => 1);

my $output;

my $peak_calling_statistics_sorted = [ sort { $a->total_length <=> $b->total_length } @$peak_calling_statistics ];

my $genome_container = Bio::EnsEMBL::Registry->get_adaptor( $species, 'core', 'GenomeContainer' );
my $mouse_ref_length = $genome_container->get_ref_length;

my $graph_display_feature_types = [
    map {
        $feature_type_adaptor->fetch_by_name($_)
    }
    (
        "H3K4me1", 
        "H3K4me2", 
        "H3K4me3", 
        "H3K9ac", 
        "H3K9me3",
        "H3K27ac", 
        "H3K27me3", 
        "H3K36me3", 
        "DNase1",
        "CTCF"
    )
];

my $graph_display_epigenomes = [
    map {
        $epigenome_adaptor->fetch_by_production_name($_)
    }
    qw(
        brain_e14_5d
        CH12_LX
        embryonic_facial_prominence_embryonic_10_5
        ES_Bruce4_embryonic__
        forebrain_embryonic_10_5
        forebrain_postnatal_0
        heart_a8w
        heart_embryonic_10_5
        heart_postnatal_0
        hindbrain_embryonic_10_5
        hindbrain_postnatal_0
        intestine_embryonic_14_5
        intestine_postnatal_0
        kidney_a8w
        kidney_embryonic_14_5
        kidney_postnatal_0
        limb_embryonic_10_5
        liver_a8w
        liver_embryonic_11_5
        liver_postnatal_0
        lung_embryonic_14_5
        lung_postnatal_0
        MEL
        MEL_cell_line____
        midbrain_embryonic_10_5
        midbrain_postnatal_0
        neural_tube_embryonic_11_5
        spleen_a8w
        stomach_embryonic_14_5
        stomach_postnatal_0
        thymus_a8w
    ),
];

my $output_file = "$output_directory/peak_calling_report.html";

$tt->process(
    $description_template, 
    {
        peak_calling_statistics => $peak_calling_statistics_sorted,
        peak_calling_adaptor    => $peak_calling_adaptor,
        
        length_to_percent => sub {
            my $length = shift;
            return $length / $mouse_ref_length;
        },
        
        round_percent => sub {
            my $number = shift;
            return sprintf("%.2f", $number);
        },

        feature_types => $graph_display_feature_types,
        epigenomes    => $graph_display_epigenomes,
    },
    "$output_directory/peak_calling_report.html"
)
    || die $tt->error;

$logger->info("Report written to $output_file\n");
$logger->finish_log;


