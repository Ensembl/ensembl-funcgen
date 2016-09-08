#!/usr/bin/env perl

use strict;
use File::Spec;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::DBSQL::DBConnection;
use Getopt::Long;

my $registry;
my $species;
my $ontology_database_url;
my $output_dir;
my $epigenome_name;
my $only_summary;

GetOptions (
   'registry=s'       => \$registry,
   'species=s'        => \$species,
   'output_dir=s'     => \$output_dir,
   'epigenome_name=s' => \$epigenome_name,
   'only_summary!'    => \$only_summary,
);

Bio::EnsEMBL::Registry->load_all($registry);

use Bio::EnsEMBL::Utils::Logger;
my $logger = Bio::EnsEMBL::Utils::Logger->new();

my $ontology_term_adaptor      = Bio::EnsEMBL::Registry->get_adaptor( 'Multi', 'Ontology', 'OntologyTerm' );
my $regulatory_feature_adaptor = Bio::EnsEMBL::Registry->get_adaptor( $species, 'Funcgen', 'RegulatoryFeature' );
my $epigenome_adaptor          = Bio::EnsEMBL::Registry->get_adaptor( $species, 'Funcgen', 'Epigenome' );

my $output_file;
my $serialiser_callback;

my $exported_something = 1;

my $parameters = {
  logger                     => $logger,
  regulatory_feature_adaptor => $regulatory_feature_adaptor, 
  epigenome_adaptor          => $epigenome_adaptor, 
  epigenome_name             => $epigenome_name,
  exported_something         => \$exported_something
};

if ($only_summary) {
  ($output_file, $serialiser_callback) = create_file_and_serialiser_for_summary($parameters);
} else {
  ($output_file, $serialiser_callback) = create_file_and_serialiser($parameters);
}

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

my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'Funcgen' );

my $helper = Bio::EnsEMBL::Utils::SqlHelper->new(
  -DB_CONNECTION => $funcgen_adaptor->dbc
);

my $last_id = 0;
my $batch_size = 10000;

while ($exported_something) {

  $exported_something = undef;

  $helper->execute_no_return(
    -SQL      => 'select regulatory_feature_id from regulatory_feature join regulatory_build using (regulatory_build_id) where is_current=1 and regulatory_feature_id > ? order by regulatory_feature_id limit ?',
    -PARAMS => [ $last_id, $batch_size ],
    -CALLBACK => $serialiser_callback,
  );
}

$logger->info("Export done.\n");
$logger->info("Gzipping $output_file\n");

use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd );
run_system_cmd("gzip $output_file");

sub create_file_and_serialiser_for_summary {

  my $param = shift;

  my $logger                     = $param->{logger};
  my $regulatory_feature_adaptor = $param->{regulatory_feature_adaptor};
  my $exported_something         = $param->{exported_something};

  my $output_file = File::Spec->catfile(
    $output_dir,
    'RegulatoryFeatures.gff'
  );
  $logger->info("The features will be written to " . $output_file ."\n");

  my $regulatory_activity_serialiser = sub {
    my @row  = @{ shift @_ };
    my $dbid = $row[0];
    my $regulatory_feature  = $regulatory_feature_adaptor->fetch_by_dbID($dbid);

    $serializer->print_feature($regulatory_feature);

    # This prevents memory leaks.
    undef %$regulatory_feature;

    $$exported_something = 1;
    $last_id = $dbid;
    return;
  };
  return ($output_file, $regulatory_activity_serialiser);
}

sub create_file_and_serialiser {

  my $param = shift;

  my $logger                     = $param->{logger};
  my $regulatory_feature_adaptor = $param->{regulatory_feature_adaptor};
  my $epigenome_adaptor          = $param->{epigenome_adaptor};
  my $epigenome_name             = $param->{epigenome_name};
  my $exported_something         = $param->{exported_something};

  my $epigenome = $epigenome_adaptor->fetch_by_name($epigenome_name);

  if (! defined $epigenome) {
    die("Can't find epigenome with name " . $epigenome_name);
  }
  $logger->info("Exporting regulatory features for " . $epigenome->display_label ."\n");

  my $output_file = File::Spec->catfile(
    $output_dir,
    $epigenome->production_name . '.gff'
  );
  $logger->info("The features will be written to " . $output_file ."\n");

  my $regulatory_activity_serialiser = sub {
    my @row  = @{ shift @_ };
    my $dbid = $row[0];
    my $regulatory_feature  = $regulatory_feature_adaptor->fetch_by_dbID($dbid);
    my $regulatory_activity = $regulatory_feature->regulatory_activity_for_epigenome($epigenome);

    if ($regulatory_activity) {
      $serializer->print_feature($regulatory_activity);
    }
    # This prevents memory leaks.
    undef %$regulatory_feature;

    $$exported_something = 1;
    $last_id = $dbid;
    return;
  };
  return ($output_file, $regulatory_activity_serialiser);
}

# sub create_filename_from_epigenome {
#   my $epigenome = shift;
#   return
#     File::Spec->catfile(
#       &gff_output_directory,
#       $epigenome->production_name
#       . '.gff'
#     )
# }
# 
# sub gff_output_directory {
#   return 'RegulatoryFeatureActivity'
# }
