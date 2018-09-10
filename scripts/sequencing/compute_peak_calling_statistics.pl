#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=head1 NAME

  A script to compute the number of bases covered by regulatory features both in total and by feature type

bsub -q production-rh7 -M8000  -R"select[mem>8000]  rusage[mem=8000]" perl scripts/sequencing/compute_peak_calling_statistics.pl \
    --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb7.pm \
    --species homo_sapiens \
    --feature_type DNase1 \
    --tempdir /hps/nobackup/production/sds-flicek-bp/blueprint_fastq_files/mnuhn_regbuild_pipeline/peak_stats

bsub -q production-rh7 -M8000  -R"select[mem>8000]  rusage[mem=8000]" perl scripts/sequencing/compute_peak_calling_statistics.pl \
    --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb7.pm \
    --species homo_sapiens \
    --feature_type CTCF \
    --tempdir /hps/nobackup/production/sds-flicek-bp/blueprint_fastq_files/mnuhn_regbuild_pipeline/peak_stats

bsub -q production-rh7 -M8000  -R"select[mem>8000]  rusage[mem=8000]" perl scripts/sequencing/compute_peak_calling_statistics.pl \
    --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb7.pm \
    --species homo_sapiens \
    --feature_type H3K4me1 \
    --tempdir /hps/nobackup/production/sds-flicek-bp/blueprint_fastq_files/mnuhn_regbuild_pipeline/peak_stats

bsub -q production-rh7 -M8000  -R"select[mem>8000]  rusage[mem=8000]" perl scripts/sequencing/compute_peak_calling_statistics.pl \
    --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb7.pm \
    --species homo_sapiens \
    --feature_type H3K4me2 \
    --tempdir /hps/nobackup/production/sds-flicek-bp/blueprint_fastq_files/mnuhn_regbuild_pipeline/peak_stats

bsub -q production-rh7 -M8000  -R"select[mem>8000]  rusage[mem=8000]" perl scripts/sequencing/compute_peak_calling_statistics.pl \
    --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb7.pm \
    --species homo_sapiens \
    --feature_type H3K4me3 \
    --tempdir /hps/nobackup/production/sds-flicek-bp/blueprint_fastq_files/mnuhn_regbuild_pipeline/peak_stats

bsub -q production-rh7 -M8000  -R"select[mem>8000]  rusage[mem=8000]" perl scripts/sequencing/compute_peak_calling_statistics.pl \
    --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb7.pm \
    --species homo_sapiens \
    --feature_type H3K9ac \
    --tempdir /hps/nobackup/production/sds-flicek-bp/blueprint_fastq_files/mnuhn_regbuild_pipeline/peak_stats

bsub -q production-rh7 -M8000  -R"select[mem>8000]  rusage[mem=8000]" perl scripts/sequencing/compute_peak_calling_statistics.pl \
    --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb7.pm \
    --species homo_sapiens \
    --feature_type H3K9me3 \
    --tempdir /hps/nobackup/production/sds-flicek-bp/blueprint_fastq_files/mnuhn_regbuild_pipeline/peak_stats

bsub -q production-rh7 -M8000  -R"select[mem>8000]  rusage[mem=8000]" perl scripts/sequencing/compute_peak_calling_statistics.pl \
    --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb7.pm \
    --species homo_sapiens \
    --feature_type H3K27ac \
    --tempdir /hps/nobackup/production/sds-flicek-bp/blueprint_fastq_files/mnuhn_regbuild_pipeline/peak_stats

bsub -q production-rh7 -M8000  -R"select[mem>8000]  rusage[mem=8000]" perl scripts/sequencing/compute_peak_calling_statistics.pl \
    --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb7.pm \
    --species homo_sapiens \
    --feature_type H3K27me3 \
    --tempdir /hps/nobackup/production/sds-flicek-bp/blueprint_fastq_files/mnuhn_regbuild_pipeline/peak_stats

bsub -q production-rh7 -M8000  -R"select[mem>8000]  rusage[mem=8000]" perl scripts/sequencing/compute_peak_calling_statistics.pl \
    --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb7.pm \
    --species homo_sapiens \
    --feature_type H3K36me3 \
    --tempdir /hps/nobackup/production/sds-flicek-bp/blueprint_fastq_files/mnuhn_regbuild_pipeline/peak_stats





=cut

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use Bio::EnsEMBL::Utils::Logger;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd run_backtick_cmd );

my %options;
GetOptions (
    \%options,
    "species|s=s",
    "registry|r=s",
    "tempdir|t=s",
    "feature_type|f=s",
 );

use Hash::Util qw( lock_keys );
lock_keys( %options );

my $species  = $options{'species'};
my $registry = $options{'registry'};
my $tempdir  = $options{'tempdir'};
my $feature_type_name = $options{'feature_type'};

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

$logger->info("species  = $species  \n");
$logger->info("registry = $registry \n");
$logger->info("tempdir  = $tempdir  \n");

Bio::EnsEMBL::Registry->load_all($registry);
my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');

my $peak_calling_statistic_adaptor = $funcgen_adaptor->get_PeakCallingStatisticAdaptor;

my $genome = Bio::EnsEMBL::Registry->get_adaptor( $species, "core", "GenomeContainer" );
my $genome_length_in_bp = $genome->get_ref_length;

$logger->info("Number of bp in $species is $genome_length_in_bp\n");

my $feature_type_adaptor = $funcgen_adaptor->get_FeatureTypeAdaptor;
my $peak_calling_adaptor = $funcgen_adaptor->get_PeakCallingAdaptor;
my $peak_adaptor         = $funcgen_adaptor->get_PeakAdaptor;

#my $feature_type_name = 'DNase1';
my $feature_type  = $feature_type_adaptor->fetch_by_name($feature_type_name);
my $peak_callings = $peak_calling_adaptor->fetch_all_by_FeatureType($feature_type);

$logger->info("Got " . @$peak_callings . " peak_callings for $feature_type_name.\n");
$logger->info("In the following epigenomes:\n");

foreach my $peak_calling (@$peak_callings) {

    my $epigenome = $peak_calling->fetch_Epigenome;
    $logger->info("     " . $epigenome->display_label . "\n");

}

my $export_directory = $tempdir . '/' . $feature_type_name;
my $all_peaks_file   = $tempdir . '/' . $feature_type_name . '.bed';

use File::Path qw( make_path );
make_path($export_directory);

$logger->info("Exporting peaks to $export_directory:\n");

open my $all_bed_fh, '>', $all_peaks_file or die("Couldn't open $all_peaks_file for writing!");

foreach my $peak_calling (@$peak_callings) {

    my $epigenome        = $peak_calling->fetch_Epigenome;
    my $export_file      = $export_directory . '/' . $epigenome->production_name . '.bed';

    $logger->info("     Writing $feature_type_name peaks for " . $epigenome->display_label . " to $export_file\n");
    
    open my $bed_fh, '>', $export_file or die("Couldn't open $export_file for writing!");
    $peak_adaptor->_bulk_export_to_bed_by_PeakCalling($peak_calling, $bed_fh, $all_bed_fh);
    $bed_fh->close;
    
    my $sort_cmd = "bedSort $export_file $export_file";
    $logger->info("     Sorting:\n");
    $logger->info("     $sort_cmd\n");
    
    run_system_cmd($sort_cmd);
}

$all_bed_fh->close;

my $sort_cmd = "bedSort $all_peaks_file $all_peaks_file";
$logger->info("     Sorting $all_peaks_file:\n");
$logger->info("     $sort_cmd\n");
run_system_cmd($sort_cmd);

$logger->info("Computing overlap\n");

my $count_cmd = "wiggletools AUC unit $all_peaks_file";
$logger->info("     Counting bases that features in $all_peaks_file overlap\n");
$logger->info("     $count_cmd\n");
my $number_of_bases = run_backtick_cmd($count_cmd);

my $coverage_percent = 100 * $number_of_bases / $genome_length_in_bp;

$logger->info("     The number of bases overlapped by $feature_type_name features is $number_of_bases\n");
$logger->info("     As a percentage of the genome this is $coverage_percent\n");

$logger->info("Storing statistics\n");

use Bio::EnsEMBL::Funcgen::PeakCallingStatistic;

$peak_calling_statistic_adaptor->store(
    Bio::EnsEMBL::Funcgen::PeakCallingStatistic->new(
        -feature_type_id => $feature_type->dbID,
        -statistic       => 'coverage_bp',
        -value           => $number_of_bases,
    )
);

$peak_calling_statistic_adaptor->store(
    Bio::EnsEMBL::Funcgen::PeakCallingStatistic->new(
        -feature_type_id => $feature_type->dbID,
        -statistic       => 'coverage_percent',
        -value           => $coverage_percent,
    )
);

$logger->info("All done.\n");
$logger->finish_log;
exit(0);

