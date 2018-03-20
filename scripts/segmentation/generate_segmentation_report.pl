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

perl scripts/segmentation/generate_segmentation_report.pl \
    --species          mus_musculus \
    --registry         /homes/mnuhn/work_dir_ersa/lib/ensembl-funcgen/registry.pm \
    --output_directory ./segmentation_report

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

my $segmentation_state_emission_adaptor = $mouse_funcgen_dba->get_SegmentationStateEmissionAdaptor;

my $segmentation_state_emissions = $segmentation_state_emission_adaptor->fetch_all;

my $file = __FILE__;
use File::Basename qw( dirname basename );
my $description_template = dirname($file) . '/../../templates/segmentation/report.html';

if (! -e $description_template) {
    die("Can't find $description_template");
}

use File::Path qw( make_path );
make_path( $output_directory );

use Template;
my $tt = Template->new( ABSOLUTE => 1, RELATIVE => 1 );

my $output;

my $segmentation_state_emissions_sorted = [ sort { $a->state <=> $b->state } @$segmentation_state_emissions ];

my $output_file = "$output_directory/segmentation_report.html";

$tt->process(
    $description_template, 
    {
        segmentation_state_emissions => $segmentation_state_emissions_sorted,
        dbc => $mouse_funcgen_dba->dbc,
        time => sub {
          return "" . localtime
        },
        round_num => sub {
            my $number = shift;
            return sprintf("%.2f", $number);
        },

    },
    $output_file
)
    || die $tt->error;

$logger->info("Report written to $output_file\n");
$logger->finish_log;


