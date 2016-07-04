#!/usr/bin/env perl

use strict;
use Bio::EnsEMBL::Registry;
use File::Spec;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::DBSQL::DBConnection;
use Getopt::Long;

my $regulation_database_url;
my $ontology_database_url;
my $ftp_base_dir;

GetOptions (
   'regulation_database_url=s' => \$regulation_database_url,
   'ontology_database_url=s'   => \$ontology_database_url,
   'ftp_base_dir=s'            => \$ftp_base_dir,
);

my $ontology_dbc   = Bio::EnsEMBL::Hive::DBSQL::DBConnection->new(-url => $ontology_database_url);
my $regulation_dbc = Bio::EnsEMBL::Hive::DBSQL::DBConnection->new(-url => $regulation_database_url);

use Bio::EnsEMBL::Utils::Logger;
my $logger = Bio::EnsEMBL::Utils::Logger->new();

use Bio::EnsEMBL::DBSQL::OntologyTermAdaptor;
my $ontology_term_adaptor = Bio::EnsEMBL::DBSQL::OntologyTermAdaptor->new(
  Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    -dbconn  => $ontology_dbc,
  )
);

my $funcgen_db_adaptor = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
    -dbconn  => $regulation_dbc,
);

my $output_file = File::Spec->catfile(
  $ftp_base_dir, 'AnnotatedFeatures.gff'
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

my $helper = Bio::EnsEMBL::Utils::SqlHelper->new(
  -DB_CONNECTION => $regulation_dbc
);

my $number_of_annotated_features = $helper->execute_simple(
  -SQL      => 'select count(annotated_feature_id) from annotated_feature',
)->[0];

$logger->info("There are " . $number_of_annotated_features ." annotated features\n");

my $annotated_feature_adaptor = $funcgen_db_adaptor->get_AnnotatedFeatureAdaptor;
my $progressbar_id = $logger->init_progress($number_of_annotated_features, 100);
my $i=0;

my $last_id = 0;
my $exported_something = 1;
my $batch_size = 100000;

while ($exported_something) {

  $exported_something = undef;

  $helper->execute_no_return(
    -SQL      => 'select annotated_feature_id from annotated_feature where annotated_feature_id > ? order by annotated_feature_id limit ?',
    -PARAMS => [ $last_id, $batch_size ],
    -CALLBACK => sub {
      my @row = @{ shift @_ };
      my $annotated_feature_id = $row[0];
      my $annotated_feature = $annotated_feature_adaptor->fetch_by_dbID($annotated_feature_id);

      $serializer->print_feature($annotated_feature);

      # This prevents memory leaks.
      undef %$annotated_feature;

      $exported_something = 1;
      $last_id = $annotated_feature_id;
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
