#!/usr/bin/env perl

use strict;
use Bio::EnsEMBL::Registry;
use File::Spec;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::DBSQL::DBConnection;
use Getopt::Long;

my $registry;
my $species;
my $output_file;
my $min_id;
my $max_id;

GetOptions (
   'registry=s'                => \$registry,
   'species=s'                 => \$species,
   'output_file=s'             => \$output_file,
   'min_id=s'                  => \$min_id,
   'max_id=s'                  => \$max_id,
);

Bio::EnsEMBL::Registry->load_all($registry);

use Bio::EnsEMBL::Utils::Logger;
my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->info("The features will be written to " . $output_file ."\n");

use File::Basename;
my $ftp_dir = dirname($output_file);

use File::Path qw(make_path);
make_path($ftp_dir);

use IO::File;
my $output_fh = IO::File->new(">$output_file");

my $ontology_term_adaptor = Bio::EnsEMBL::Registry->get_adaptor( 'Multi', 'Ontology', 'OntologyTerm' );

use Bio::EnsEMBL::Utils::IO::GFFSerializer;
my $serializer = Bio::EnsEMBL::Utils::IO::GFFSerializer->new(
  $ontology_term_adaptor,
  $output_fh
);

my $annotated_feature_adaptor = Bio::EnsEMBL::Registry->get_adaptor( $species, 'Funcgen', 'AnnotatedFeature' );
my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'Funcgen' );

my $helper = Bio::EnsEMBL::Utils::SqlHelper->new(
  -DB_CONNECTION => $funcgen_adaptor->dbc
);

my $batch_constraint = qq(annotated_feature_id<=$max_id and annotated_feature_id>=$min_id);

my $number_of_annotated_features = $helper->execute_simple(
  -SQL      => "select count(annotated_feature_id) from annotated_feature where $batch_constraint",
)->[0];

$logger->info("About to export " . $number_of_annotated_features ." annotated features\n");

my $progressbar_id = $logger->init_progress($number_of_annotated_features, 100);
my $i=0;

my $last_id = 0;
my $exported_something = 1;
my $batch_size = 1000;

while ($exported_something) {

  $exported_something = undef;

  $helper->execute_no_return(
    -SQL      => "select annotated_feature_id from annotated_feature where annotated_feature_id > ? and $batch_constraint order by annotated_feature_id limit ?",
    -PARAMS => [ $last_id, $batch_size ],
    -CALLBACK => sub {
      my @row = @{ shift @_ };
      my $annotated_feature_id = $row[0];
      my $annotated_feature = $annotated_feature_adaptor->fetch_by_dbID($annotated_feature_id);

      eval {
	$serializer->print_feature($annotated_feature);
      };
      if ($@) {
	use Carp;
	confess(
	  "Unable to serialise feature! dbid:"
	  . $annotated_feature->dbID
	);
      }

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
# $logger->info("Gzipping $output_file\n");
# 
# use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd );
# run_system_cmd("gzip $output_file");
