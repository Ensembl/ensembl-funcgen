#!/usr/bin/env perl

use strict;
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use Carp;
use Bio::EnsEMBL::Utils::Logger;

=head1 


delete from data_file where table_name = "segmentation_file";
truncate segmentation_file;

select * from data_file where table_name = "segmentation_file"

rm -rf /gpfs/nobackup/ensembl/mnuhn/mnuhn/regulatory_build_pipeline_run3/dbfiles/mus_musculus/GRCm38/funcgen/segmentation_file

perl scripts/segmentation/register_segmentation_files.pl  \
    --species                          mus_musculus \
    --registry                         /homes/mnuhn/work_dir_dev_break/lib/ensembl-funcgen/registry.with_previous_version.pm \
    --projected_segmentation_directory /gpfs/nobackup/ensembl/mnuhn/mnuhn/regulatory_build_pipeline_run3/temp_dir/regulatory_build/mus_musculus/GRCm38/projected_segmentations/ \
    --db_file_species_assembly_dir     /gpfs/nobackup/ensembl/mnuhn/mnuhn/regulatory_build_pipeline_run3/dbfiles/mus_musculus/GRCm38/funcgen/segmentation_file/093 \
    --db_file_relative_dir             /funcgen/segmentation_file/093

=cut

use strict;
use Getopt::Long;

my $species;
my $registry;
my $output_directory;
my $projected_segmentation_directory;
my $db_file_species_assembly_dir;
my $db_file_relative_dir;

GetOptions (
   'species=s'                          => \$species,
   'registry=s'                         => \$registry,
   'projected_segmentation_directory=s' => \$projected_segmentation_directory,
   'db_file_species_assembly_dir=s'     => \$db_file_species_assembly_dir,
   'db_file_relative_dir=s'             => \$db_file_relative_dir,
);

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

$logger->info("registry                          = " . $registry         . "\n");
$logger->info("species                           = " . $species          . "\n");
$logger->info("projected_segmentation_directory  = " . $projected_segmentation_directory . "\n");
$logger->info("db_file_species_assembly_dir      = " . $db_file_species_assembly_dir . "\n");
$logger->info("db_file_relative_dir              = " . $db_file_relative_dir . "\n");

use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->load_all($registry);

my $funcgen_dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');

my $regulatory_build_adaptor  = $funcgen_dba->get_RegulatoryBuildAdaptor;
my $segmentation_file_adaptor = $funcgen_dba->get_SegmentationFileAdaptor;
my $analysis_adaptor          = $funcgen_dba->get_AnalysisAdaptor;

my $current_regulatory_build = $regulatory_build_adaptor->fetch_current_regulatory_build;
my $all_epigenomes = $current_regulatory_build->get_all_Epigenomes;

use File::Path qw( make_path );
make_path($db_file_species_assembly_dir);

my $segmentation_analysis = $analysis_adaptor->fetch_by_dbID(1);

foreach my $current_epigenome (@$all_epigenomes) {

  my $big_bed_file_basename = $current_epigenome->production_name . '.bb';
  my $big_bed_file = $projected_segmentation_directory . '/' . $big_bed_file_basename;
  
  if (! -e $big_bed_file) {
    die("$big_bed_file doesn't exist!");
  }

  $logger->info("Computing checksum\n");
  
  use Digest::MD5 qw(md5 md5_hex md5_base64);
  open (my $fh, '<', $big_bed_file) or die "Can't open '$big_bed_file': $!";
  binmode ($fh);
  my $md5sum = Digest::MD5->new->addfile($fh)->hexdigest;
  $fh->close;
  
  $logger->info("$md5sum\n");

  my $destination = $db_file_species_assembly_dir . '/' . $big_bed_file_basename;

  $logger->info("Copying $big_bed_file to $destination\n");
  
  my $registered_name = $db_file_relative_dir . '/' . $big_bed_file_basename;
  
#   my $cmd = qq(cp $big_bed_file $destination);
  use File::Copy;
  copy($big_bed_file, $destination) || die("Couldn't copy $big_bed_file to $destination!");
  
  $logger->info("Registering $registered_name\n");
  
  use Bio::EnsEMBL::Funcgen::SegmentationFile;
  
  my $segmentation_file = Bio::EnsEMBL::Funcgen::SegmentationFile->new(
    -name             => 'Segmentation of ' . $current_epigenome->name,
    -analysis         => $segmentation_analysis,
    -epigenome        => $current_epigenome,
    -regulatory_build => $current_regulatory_build,
    -file             => $registered_name,
    -file_type        => 'BIGBED',
    -md5sum           => $md5sum,
  );
  $segmentation_file_adaptor->store($segmentation_file);
}

$logger->finish_log;


