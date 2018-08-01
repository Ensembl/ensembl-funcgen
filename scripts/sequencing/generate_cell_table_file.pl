#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Logger;

my $registry;
my $species;
my $file;
my $class;
my $superclass;

=head1

generate_cell_table_file.pl \
    --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb4.pm \
    --species homo_sapiens \
    --class no_ctcf \
    --superclass blueprint \
    --file deleteme.txt

=cut

my %config_hash = (
  'registry'   => \$registry,
  'species'    => \$species,
  'cell_table_file'       => \$file,
  'class'      => \$class,
  'superclass' => \$superclass,
);

my $result = GetOptions(
  \%config_hash,
  'registry=s',
  'species=s',
  'cell_table_file=s',
  'class=s',
  'superclass=s',
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
    epigenome,
    feature_type,
    signal_bam_path,
    control_bam_path 
  from 
    segmentation_cell_tables 
  where 
    superclass = "$superclass" 
    and class  = "$class"
SQL
;

use File::Basename qw( dirname fileparse );
use File::Path qw(make_path remove_tree);

my $output_dir = dirname($file);
make_path($output_dir);

open my $fh, '>', $file or die("Can't open $class");

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
$logger->info("Written to $file\n");
$logger->finish_log;
