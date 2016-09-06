#!/usr/bin/env perl

use strict;
use Bio::EnsEMBL::Registry;
use File::Spec;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::DBSQL::DBConnection;
use Getopt::Long;

my $registry;
my $species;
my $ftp_base_dir;

GetOptions (
   'registry=s'                => \$registry,
   'species=s'                 => \$species,
   'ftp_base_dir=s'            => \$ftp_base_dir,
);

Bio::EnsEMBL::Registry->load_all($registry);

use Bio::EnsEMBL::Utils::Logger;
my $logger = Bio::EnsEMBL::Utils::Logger->new();

my $ontology_term_adaptor = Bio::EnsEMBL::Registry->get_adaptor( 'Multi', 'Ontology', 'OntologyTerm' );

my $output_file = File::Spec->catfile(
  $ftp_base_dir, 'MotifFeatures.gff'
);
$logger->info("The features will be written to " . $output_file ."\n");

use File::Basename;
my $ftp_dir = dirname($output_file);

use File::Path qw(make_path);
make_path($ftp_dir);

use IO::File;
my $output_fh = IO::File->new(">$output_file");

use Bio::EnsEMBL::Utils::IO::GFFSerializer;
my $serializer = Bio::EnsEMBL::Utils::IO::GFFSerializer->new(
  $ontology_term_adaptor,
  $output_fh
);

my $motif_feature_adaptor = Bio::EnsEMBL::Registry->get_adaptor( $species, 'Funcgen', 'MotifFeature' );
my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'Funcgen' );

my $helper = Bio::EnsEMBL::Utils::SqlHelper->new(
  -DB_CONNECTION => $funcgen_adaptor->dbc
);

my $number_of_motif_features = $helper->execute_simple(
  -SQL      => 'select count(motif_feature_id) from motif_feature',
)->[0];

$logger->info("There are " . $number_of_motif_features ." motif features\n");

my $progressbar_id = $logger->init_progress($number_of_motif_features, 100);
my $i=0;

my $last_id = 0;
my $exported_something = 1;
my $batch_size = 100000;

while ($exported_something) {

  $exported_something = undef;

  $helper->execute_no_return(
    -SQL      => 'select motif_feature_id from motif_feature where motif_feature_id > ? order by motif_feature_id limit ?',
    -PARAMS => [ $last_id, $batch_size ],
    -CALLBACK => sub {
      my @row = @{ shift @_ };
      my $motif_feature_id = $row[0];
      my $motif_feature = $motif_feature_adaptor->fetch_by_dbID($motif_feature_id);

      $serializer->print_feature($motif_feature);

      # This prevents memory leaks.
      undef %$motif_feature;

      $exported_something = 1;
      $last_id = $motif_feature_id;
      return;
    },
  );
  $i+=$batch_size;
  $logger->log_progressbar($progressbar_id, $i);

}
$logger->info("Export done.\n");
$logger->info("Gzipping $output_file\n");

use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd );
run_system_cmd("gzip $output_file");
