#!/usr/bin/env perl

use strict;
use JSON;
use Bio::EnsEMBL::Registry;
use File::Spec;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::DBSQL::DBConnection;
use Getopt::Long;

# export_qc_proportion_of_reads_in_peaks.pl --registry /nfs/users/nfs_m/mn1/work_dir_ftp/lib/ensembl-funcgen/registry.pm --species homo_sapiens --output_file testun/qc_proportion_of_reads_in_peaks.json 
# export_qc_proportion_of_reads_in_peaks.pl --registry /nfs/users/nfs_m/mn1/work_dir_ftp/lib/ensembl-funcgen/registry.pm --species homo_sapiens | less

=head1

export_qc_proportion_of_reads_in_peaks.pl \
  --output_file /hps/nobackup/production/ensembl/mnuhn/chip_seq_analysis/ftp/homo_sapiens/QualityChecks/homo_sapiens.GRCh38.proportion_of_reads_in_peaks.quality_check.20171220.json  \
  --registry /homes/mnuhn/work_dir_ersa/lib/ensembl-funcgen/registry.pm  \
  --species homo_sapiens

=cut

my $registry;
my $species;
my $output_file;

GetOptions (
   'registry=s'    => \$registry,
   'species=s'     => \$species,
   'output_file=s' => \$output_file,
);

Bio::EnsEMBL::Registry->load_all($registry);

use Bio::EnsEMBL::Utils::Logger;
my $logger = Bio::EnsEMBL::Utils::Logger->new();

my $funcgen_db_adaptor   = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
my $phantom_peak_adaptor = $funcgen_db_adaptor->get_PhantomPeakAdaptor;

my $json = JSON->new->utf8;
$json->pretty(1);
$json->canonical(1);

my $output_fh;
if ($output_file) {
  $logger->info("The features will be written to " . $output_file ."\n");

  use File::Basename;
  my $ftp_dir = dirname($output_file);

  use File::Path qw(make_path);
  make_path($ftp_dir);

  use IO::File;
  $output_fh = IO::File->new(">$output_file");
} else {
  $output_fh = *STDOUT;
}

my $all_phantom_peak_objects = $phantom_peak_adaptor->fetch_all;
#my $all_phantom_peak_objects = $phantom_peak_adaptor->fetch_all('LIMIT 5');

my $phantom_peak_summaries = [ map { $_->summary_as_hash } @$all_phantom_peak_objects ];

$output_fh->print(
  $json->encode($phantom_peak_summaries)
);
