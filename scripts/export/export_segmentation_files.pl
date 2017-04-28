#!/usr/bin/env perl

use strict;
use Bio::EnsEMBL::Registry;
use File::Spec;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::DBSQL::DBConnection;
use Getopt::Long;
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::Utils::ExportUtils qw(
  assert_source_files_exist
  assert_destination_file_names_uniqe
);

=head1

export_segmentation_files.pl \
  --registry /nfs/users/nfs_m/mn1/work_dir_ftp/lib/ensembl-funcgen/registry.pm  \
  --species homo_sapiens       \
  --assembly GrCh38            \
  --data_freeze_date 20161007  \
  --dbfile_registry_path /nfs/ensnfs-webdev/staging/homo_sapiens/GRCh38  \
  --dbfile_registry_path /lustre/scratch109/ensembl/funcgen/mn1/ersa/mn1_homo_sapiens_funcgen_86_38/output/mn1_homo_sapiens_funcgen_86_38/funcgen  \
  --dbfile_registry_path /lustre/scratch109/ensembl/funcgen/mn1/ersa/mn1_dev3_homo_sapiens_funcgen_85_38/output/mn1_dev3_homo_sapiens_funcgen_85_38/funcgen  \
  --destination_root_path /lustre/scratch109/ensembl/funcgen/mn1/ftpsite10/homo_sapiens/Segmentation

=cut

my $registry;
my $species;
my $assembly;
my $data_freeze_date;
my $destination_root_path;
my @dbfile_registry_path;
my $die_if_source_files_missing = 1;

GetOptions (
   'registry=s'              => \$registry,
   'species=s'               => \$species,
   'assembly=s'              => \$assembly,
   'data_freeze_date=s'      => \$data_freeze_date,
   'destination_root_path=s' => \$destination_root_path,
   'die_if_source_files_missing=s'          => \$die_if_source_files_missing,
   'dbfile_registry_path=s' => \@dbfile_registry_path,
);

Bio::EnsEMBL::Registry->load_all($registry);

use Bio::EnsEMBL::Utils::Logger;
my $logger = Bio::EnsEMBL::Utils::Logger->new();

my $segmentation_file_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'SegmentationFile');

my $all_segmentation_files = $segmentation_file_adaptor->fetch_all;

$logger->info("There are " . @$all_segmentation_files . " segmentation files in the database.\n");

my %source_file_to_destination_file_map;
SEGMENTATION_FILE:
foreach my $current_segmentation_file (@$all_segmentation_files) {

  my $epigenome_production_name        = $current_segmentation_file->get_Epigenome->production_name;
  my $analysis_logic_name              = $current_segmentation_file->get_Analysis->logic_name;
  my $is_from_current_regulatory_build = $current_segmentation_file->get_RegulatoryBuild->is_current;
  
  # Only include segmentation files from the current regulatory build
  next (SEGMENTATION_FILE) if (!$is_from_current_regulatory_build );
  
  # Find the right base path among the ones provided, if none can be found, leave empty.
  #
  my $this_files_dbfile_registry_path;
  POSSIBLE_ROOT_BASE_PATH:
  foreach my $current_dbfile_registry_path (@dbfile_registry_path) {
    my $source_file_candidate = $current_dbfile_registry_path . '/' . $current_segmentation_file->file;
    if (-e $source_file_candidate) {
      $this_files_dbfile_registry_path = $current_dbfile_registry_path;
      last POSSIBLE_ROOT_BASE_PATH;
    }
  }
  my $source_file = $this_files_dbfile_registry_path . '/' . $current_segmentation_file->file;
  
  my $destination_file_name = join '.', (
    $species,
    $assembly,
    $epigenome_production_name,
    $analysis_logic_name,
    $data_freeze_date,
    'bb'
  );
  my $destination_file = $destination_root_path . '/' .$destination_file_name;
  
  $source_file_to_destination_file_map{$source_file} = $destination_file;
}

my @non_existing_source_files = assert_source_files_exist(keys %source_file_to_destination_file_map);

if (@non_existing_source_files) {
  if ($die_if_source_files_missing) {
    die(
      "The following files from the database do not exist: " . Dumper(@non_existing_source_files) . "\n"
    );
  } else {
    $logger->warning(
      "The following files from the database do not exist: " . Dumper(@non_existing_source_files) . "\n"
    );
  }
}

$logger->info("Generating commands for creating the ftp site\n");

my @cmd;
push @cmd, qq(mkdir -p $destination_root_path);

foreach my $current_source_file (keys %source_file_to_destination_file_map) {

  my $destination_file = $source_file_to_destination_file_map{$current_source_file};
  push @cmd, qq(ln -s $current_source_file $destination_file);
}

$logger->info("Running commands for creating the ftp site\n");
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd );
foreach my $current_cmd (@cmd) {
  $logger->info($current_cmd . "\n");
  run_system_cmd($current_cmd, undef, 1);
}

$logger->info("All done.\n");

# use Data::Dumper;
# print Dumper(\%source_to_destination_file_name);

