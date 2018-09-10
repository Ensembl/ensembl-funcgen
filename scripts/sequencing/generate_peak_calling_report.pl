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

generate_peak_calling_report.pl \
    --species          mouse_with_regbuild \
    --registry         /homes/mnuhn/work_dir_ersa/lib/ensembl-funcgen/registry.pm \
    --output_directory ./reports/


generate_peak_calling_report.pl \
    --species          homo_sapiens \
    --registry         /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb7.pm \
    --output_directory /homes/mnuhn/public_html/regulatory_build_stats/rb_human_merged_old_and_new/homo_sapiens/

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
my $idr_adaptor          = $mouse_funcgen_dba->get_IdrAdaptor;

my $file = __FILE__;
use File::Basename qw( dirname basename );
my $description_template = dirname($file) . '/../../templates/peak_calling/report.html';

if (! -e $description_template) {
    die("Can't find $description_template");
}

use File::Path qw( make_path );
make_path( $output_directory );

use Template;
my $tt = Template->new( ABSOLUTE => 1, RELATIVE => 1);

my $output;

my $genome_container = Bio::EnsEMBL::Registry->get_adaptor( $species, 'core', 'GenomeContainer' );
my $genome_size_in_bp = $genome_container->get_ref_length;

my $graph_display_feature_types = [
    map {
        $feature_type_adaptor->fetch_by_name($_) || die ("Can't fetch $_");
    }
    (
        "CTCF",
        "DNase1",
        "H3K4me1", 
        "H3K4me2", 
        "H3K4me3", 
        "H3K9ac", 
        "H3K9me3",
        "H3K27ac", 
        "H3K27me3", 
        "H3K36me3", 
    )
];

my $graph_display_epigenomes = $epigenome_adaptor->fetch_all;

use Number::Format qw( format_number );

my $de = new Number::Format(
    -thousands_sep   => ',',
    -decimal_point   => '.',
);

my $output_file = "$output_directory/peak_calling_report.html";

my $dbc = $mouse_funcgen_dba->dbc;

my $experiment_adaptor = $mouse_funcgen_dba->get_ExperimentAdaptor;
my @signal_experiments  = $experiment_adaptor->_fetch_all_signal_experiments;
my @control_experiments = $experiment_adaptor->_fetch_all_control_experiments;

$Template::Stash::PRIVATE = undef;

$tt->process(
    $description_template, 
    {

        signal_experiments  => \@signal_experiments,
        control_experiments => \@control_experiments,

 #       peak_calling_statistics => $peak_calling_statistics_sorted,
        peak_calling_statistic_adaptor => $peak_calling_statistic_adaptor,
        peak_calling_adaptor    => $peak_calling_adaptor,
        idr_adaptor             => $idr_adaptor,
        dbc => $dbc,
        
        genome_size_in_bp => $genome_size_in_bp,
        
        length_to_percent => sub {
            my $length = shift;
            return $length * 100 / $genome_size_in_bp;
        },
        
        round_percent => sub {
            my $number = shift;
            return sprintf("%.2f", $number);
        },
        boolean_to_yes_no => sub {
            my $boolean = shift;
            if ($boolean) {
                return 'yes'
            }
            return 'no'
        },
        time => sub {
          return "" . localtime
        },
        feature_types => $graph_display_feature_types,
        epigenomes    => $graph_display_epigenomes,

        fetch_idr => sub  {
          my $experiment = shift;
          return $experiment->_fetch_Idr
        },
        bytes_to_gb => sub {
          my $size_in_bytes = shift;
          return ( 0 + $size_in_bytes ) / (1024 * 1024 * 1024)
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
        return $de->format_number($number);
      },


    },
    "$output_directory/peak_calling_report.html"
)
    || die $tt->error;

$logger->info("Report written to $output_file\n");
$logger->finish_log;


