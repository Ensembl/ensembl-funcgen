#!/usr/bin/env perl

use strict;
use JSON;
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

=head1 Examples

export_data_files.pl  \
  --destination_root_path /hps/nobackup/production/ensembl/mnuhn/otar/ftp/homo_sapiens/Alignments  \
  --file_type bam  \
  --assembly GRCh38 \
  --dbfile_registry_path /hps/nobackup/production/ensembl/mnuhn/otar/dbfiles/funcgen/homo_sapiens/GRCh38 \
  --registry /homes/mnuhn/work_dir_ersa/lib/ensembl-funcgen/registry.pm  \
  --die_if_source_files_missing 0 \
  --data_freeze_date 20170425 \
  --species homo_sapiens

export_data_files.pl  \                  --destination_root_path #ftp_base_dir#/#species#/Alignments  \                  --file_type bam  \                  --assembly #assembly# \                  --dbfile_registry_path /hps/nobackup/production/ensembl/mnuhn/otar/dbfiles/funcgen/#species#/#assembly# \                  --registry #reg_conf#  \                  --die_if_source_files_missing 0 \                  --data_freeze_date #data_freeze_date# \                  --species #species#                

export_data_files.pl  \                  --destination_root_path #ftp_base_dir#/#species#/Alignments  \                  --file_type bigwig  \                  --assembly #assembly# \                  --dbfile_registry_path /hps/nobackup/production/ensembl/mnuhn/otar/dbfiles/funcgen/#species#/#assembly# \                  --registry #reg_conf#  \                  --die_if_source_files_missing 0 \                  --data_freeze_date #data_freeze_date# \                  --species #species#                

=cut

my $registry;
my $species;
my @dbfile_registry_path;
my $destination_root_path;
my $file_type;
my $die_if_source_files_missing = 1;
my $assembly;

use Date::Format;
# Looks like: 20160928
my $data_freeze_date = time2str('%Y%m%d', time);

GetOptions (
   'registry=s'                    => \$registry,
   'species=s'                     => \$species,
   'file_type=s'                   => \$file_type,
   'assembly=s'                    => \$assembly,
   'data_freeze_date=s'            => \$data_freeze_date,
   'dbfile_registry_path=s'        => \@dbfile_registry_path,
   'die_if_source_files_missing=s' => \$die_if_source_files_missing,
   'destination_root_path=s'       => \$destination_root_path,
);

my $make_copies = 1;

my $creation_command = 'ln -s';
if ($make_copies) {
  $creation_command = 'cp';
}

use Bio::EnsEMBL::Utils::Logger;
my $logger = Bio::EnsEMBL::Utils::Logger->new();

if (! defined $assembly) {
  die("assembly (used in the file name only) has not been set!");
}
if (@dbfile_registry_path==0) {
  die("dbfile_registry_path has not been set!");
}
# foreach my $current_path (@dbfile_registry_path) {
#   if (! -d $current_path) {
#     $logger->error("$current_path is not an existing directory!");
#   }
# }
# if (! -d $dbfile_registry_path[0]) {
#   die($dbfile_registry_path[0] . " is not an existing directory!");
# }
if (! defined $file_type) {
  die("file_type has not been set!");
}
if ($file_type ne 'bam' && $file_type ne 'bigwig') {
  die("file_type must be one of 'bam' or 'bigwig'!");
}

Bio::EnsEMBL::Registry->load_all($registry);

my $funcgen_db_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');

my $db_connection = $funcgen_db_adaptor->dbc;

use Bio::EnsEMBL::Funcgen::Ftp::FetchQCRelatedData;
my $fetchQCRelatedData = Bio::EnsEMBL::Funcgen::Ftp::FetchQCRelatedData->new(
  -db_connection => $db_connection,
);

my $db_file_type;
my $destination_file_extension;

if ($file_type eq 'bam')    {
  $db_file_type = 'BAM';
  $destination_file_extension = 'bam';
}
if ($file_type eq 'bigwig') {
  $db_file_type = 'BIGWIG';
  $destination_file_extension = 'bw';
}

my $sth = $db_connection->prepare(
  qq(
select 
  signal_file.path signal_file, 
  result_set.name signal_alignment_name,
  feature_type.name signal_feature_type,
  control_file.path control_file,
  control_alignment.name control_alignment_name,
  control_feature_type.name control_feature_type,
  analysis.logic_name analysis,
  epigenome.display_label as epigenome,
  epigenome.production_name as epigenome_production_name
from 
  result_set
  join epigenome using (epigenome_id)
  join feature_type using (feature_type_id)
  join dbfile_registry signal_file    on (signal_file.table_name='result_set'    and signal_file.table_id=result_set_id    and signal_file.file_type='$db_file_type')
  join analysis on (result_set.analysis_id=analysis.analysis_id)
  join experiment the_signal using (experiment_id)
  left join experiment control on (the_signal.control_id=control.experiment_id)
  join result_set control_alignment on (control.experiment_id=control_alignment.experiment_id)
  join feature_type control_feature_type on (control_alignment.feature_type_id=control_feature_type.feature_type_id)
  left join dbfile_registry control_file    on (control_file.table_name='result_set'    and control_file.table_id=control_alignment.result_set_id    and control_file.file_type='$db_file_type')
order by
  epigenome_production_name,
  signal_feature_type;
  )
);
$sth->execute;

# my $output_fh;
# if ($output_file) {
#   $logger->info("The features will be written to " . $output_file ."\n");
# 
#   use File::Basename;
#   my $ftp_dir = dirname($output_file);
# 
#   use File::Path qw(make_path);
#   make_path($ftp_dir);
# 
#   use IO::File;
#   $output_fh = IO::File->new(">$output_file");
# } else {
#   $output_fh = *STDOUT;
# }

my %source_file_to_destination_file_map;
my @file_not_found;
while (my $hash_ref = $sth->fetchrow_hashref) {

  translate($hash_ref);

  # File name components as specified here:
  # https://github.com/FAANG/faang-metadata/blob/master/docs/faang_analysis_metadata.md#file-naming
  #
  my $sample_name             = $hash_ref->{epigenome_production_name};
  my $assay_type              = $hash_ref->{signal_alignment}->{feature_type};
  my $analysis_protocol_name  = $hash_ref->{analysis};
  
  $logger->info("Fetching data for " . $hash_ref->{signal_alignment}->{file} . "\n");
  
  # Find the right base path among the ones provided, if none can be found, leave empty.
  #
  my $this_files_species_assembly_data_file_base_path;
  my $found_file_in_path = undef;
  POSSIBLE_ROOT_BASE_PATH:
  foreach my $current_species_assembly_data_file_base_path (@dbfile_registry_path) {
    my $source_file_candidate = 
      $current_species_assembly_data_file_base_path 
#       . '/' . $species
#       . '/' . $assembly
      . '/' . $hash_ref->{signal_alignment}->{file};
      
    if (-e $source_file_candidate) {
      $this_files_species_assembly_data_file_base_path = $current_species_assembly_data_file_base_path;
      $found_file_in_path = 1;
      last POSSIBLE_ROOT_BASE_PATH;
    }
  }
  if (!$found_file_in_path) {
    push @file_not_found, $hash_ref->{signal_alignment}->{file};
  }
  my $source_file = $this_files_species_assembly_data_file_base_path . '/' . $hash_ref->{signal_alignment}->{file};
  
  my $replicate_description_string = create_replicate_description_string($hash_ref->{signal_alignment}->{name});
  my $destination_directory = join '/', $sample_name, $assay_type;
  my $destination_file = join '.', (
    $species,
    $assembly,
    $sample_name,
    $assay_type,
    $replicate_description_string,
    $analysis_protocol_name,
    $data_freeze_date,
    $destination_file_extension
  );
  $source_file_to_destination_file_map{$source_file} = "$destination_root_path/$destination_directory/$destination_file";
}

if (@file_not_found) {
  $logger->error(
    "The following files were not found:\n" . join "\n", map { '  - ' . $_ } @file_not_found
  );
  die;
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

foreach my $current_not_existing_file (@non_existing_source_files) {
  delete $source_file_to_destination_file_map{$current_not_existing_file};
}

my $number_of_destination_files = scalar values %source_file_to_destination_file_map;

$logger->info("Planning for " . $number_of_destination_files . " destination files. Checking that the names are unique\n");
assert_destination_file_names_uniqe(\%source_file_to_destination_file_map);

$logger->info("Generating commands for creating the ftp site\n");
# Generate a list of commands for creating directories and linking the files
my @cmd;
my %created_directories;
foreach my $current_source_file (keys %source_file_to_destination_file_map) {

  my $destination_file = $source_file_to_destination_file_map{$current_source_file};
  use File::Basename;
  my $destination_directory = dirname($destination_file);
  
  if (! exists $created_directories{$destination_directory}) {
    push @cmd, qq(mkdir -p $destination_directory);
    $created_directories{$destination_directory} = 1;
  }
  push @cmd, qq($creation_command $current_source_file $destination_file);
}

$logger->info("Running ".scalar @cmd." commands for creating the ftp site\n");
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd );
foreach my $current_cmd (@cmd) {
  $logger->info("$current_cmd\n");
  run_system_cmd($current_cmd, undef, 1);
}

$logger->info("All done.\n");

sub create_replicate_description_string {

  my $result_set_name = shift;
  
  my $result_set_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'ResultSet');
  my $result_set = $result_set_adaptor->fetch_by_name($result_set_name);

  my $support = $result_set->get_support;
  my @signal_support = grep { $_->is_control == 0 } @$support;
  
  use Bio::EnsEMBL::Funcgen::Utils::ERSAShared qw( create_replicate_input_subset_string );
  
  my $replicate_input_subset_string = create_replicate_input_subset_string(@signal_support);
  return $replicate_input_subset_string;
}

sub translate {
  my $hash_ref = shift;
  
  $hash_ref->{control_alignment} = {
    name         => $hash_ref->{control_alignment_name},
    file         => $hash_ref->{control_file},
    feature_type => $hash_ref->{control_feature_type},
  };
  delete $hash_ref->{control_alignment_name};
  delete $hash_ref->{control_file};
  delete $hash_ref->{control_feature_type};

  $hash_ref->{signal_alignment} = {
    name         => $hash_ref->{signal_alignment_name},
    file         => $hash_ref->{signal_file},
    feature_type => $hash_ref->{signal_feature_type},
  };
  delete $hash_ref->{signal_alignment_name};
  delete $hash_ref->{signal_file};
  delete $hash_ref->{signal_feature_type};
}
