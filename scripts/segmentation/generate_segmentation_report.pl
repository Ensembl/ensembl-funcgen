#!/usr/bin/env perl

use strict;
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use Carp;
use Bio::EnsEMBL::Utils::Logger;

=head1 

perl scripts/segmentation/generate_segmentation_report.pl \
    --species          homo_sapiens \
    --registry         /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb7.pm \
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

my $funcgen_dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');

my $segmentation_state_emission_adaptor = $funcgen_dba->get_SegmentationStateEmissionAdaptor;

my $segmentation_state_emissions = $segmentation_state_emission_adaptor->fetch_all;

my $file = __FILE__;
use File::Basename qw( dirname basename );
my $template_dir = dirname($file) . '/../../templates/segmentation';
my $description_template = $template_dir . '/report.html';

if (! -e $description_template) {
    die("Can't find $description_template");
}

use File::Path qw( make_path );
make_path( $output_directory );

use Template;
my $tt = Template->new( ABSOLUTE => 1, RELATIVE => 1, INCLUDE_PATH => $template_dir );

my $output;

my $segmentation_state_emissions_sorted = [ sort { $a->state <=> $b->state } @$segmentation_state_emissions ];

my $output_file = "$output_directory/segmentation_report.html";

my $segmentation_statistic_adaptor = $funcgen_dba->get_SegmentationStatisticAdaptor;
my $segmentation_adaptor           = $funcgen_dba->get_SegmentationAdaptor;

my @segmentation_assignments = qw(
  ctcf
  dead
  distal
  gene
  poised
  proximal
  repressed
  tss
  weak
);

use Number::Format qw( format_number );

my $de = new Number::Format(
    -thousands_sep   => ',',
    -decimal_point   => '.',
);

use Scalar::Util qw( looks_like_number );

$tt->process(
    $description_template, 
    {

        segmentation_statistic_adaptor 
          => $segmentation_statistic_adaptor,

        segmentation_adaptor 
          => $segmentation_adaptor,

        segmentation_assignments
          => \@segmentation_assignments,
          
        segmentation_state_emissions => $segmentation_state_emissions_sorted,
        dbc => $funcgen_dba->dbc,
        time => sub {
          return "" . localtime
        },
        round_num => sub {
            my $number = shift;
            return sprintf("%.2f", $number);
        },
        format_number => sub {
          my $number = shift;
          if (! defined $number) {
            return '-'
          }
          if ($number eq '') {
            return '-'
          }
          if ($number eq 'not defined') {
            return 0
          }
          
          if (! looks_like_number($number) ) {
            return $number;
          }
          #return 'foo';
          return $de->format_number($number);
        },

    },
    $output_file
)
    || die $tt->error;

$logger->info("Report written to $output_file\n");
$logger->finish_log;


