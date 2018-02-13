#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Logger;

my $registry;
my $species;
my $cell_table_file;

my %config_hash = (
  'registry'        => \$registry,
  'species'         => \$species,
  'cell_table_file' => \$cell_table_file,
);

my $result = GetOptions(
  \%config_hash,
  'registry=s',
  'species=s',
  'cell_table_file=s',
);

Bio::EnsEMBL::Registry->load_all($registry);

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

my $adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
my $dbc = $adaptor->dbc;

use Bio::EnsEMBL::Utils::SqlHelper;
my $sql_helper = Bio::EnsEMBL::Utils::SqlHelper->new(
  -DB_CONNECTION => $dbc
);

my $sql = <<SQL
  select
    epigenome.production_name as epigenome,
    feature_type.name as feature_type,
    signal_bam.path,
    control_bam.path
  from
    experiment
    join epigenome using (epigenome_id)
    join feature_type using (feature_type_id)
    join alignment signal_alignment on (
      experiment.experiment_id = signal_alignment.experiment_id 
      and signal_alignment.is_complete=1 
      and signal_alignment.has_duplicates=0
    )
    join data_file signal_bam on (signal_alignment.bam_file_id = signal_bam.data_file_id)
    join experiment control_experiment on (control_experiment.experiment_id = experiment.control_id)
    join alignment control_alignment on (
      control_experiment.experiment_id = control_alignment.experiment_id 
      and control_alignment.is_complete=1 
      and control_alignment.has_duplicates=0
    )
    join data_file control_bam on (control_alignment.bam_file_id = control_bam.data_file_id)
  where 
    feature_type.name in (
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
  order by 
    epigenome, feature_type
SQL
;

use File::Basename qw( dirname fileparse );
use File::Path qw(make_path remove_tree);

my $output_dir = dirname($cell_table_file);
make_path($output_dir);

open my $fh, '>', $cell_table_file or die("Can't open $cell_table_file");

$sql_helper->execute_no_return(
  -SQL          => $sql,
  -USE_HASHREFS => 0,
  -CALLBACK     => sub {
      my $row = shift;
      $fh->print(join "\t", @$row);
      $fh->print("\n");
      return;
    },
);

$fh->close;
$logger->info("Written to $cell_table_file\n");
$logger->finish_log;
