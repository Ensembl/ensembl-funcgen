#!/usr/bin/env perl

use strict;
use File::Spec;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::DBSQL::DBConnection;
use Getopt::Long;

my $registry;
my $species;
my $ontology_database_url;
my $output_file;
my $epigenome_name;
my $only_summary;

GetOptions (
   'registry=s'       => \$registry,
   'species=s'        => \$species,
   'output_file=s'     => \$output_file,
   'epigenome_name=s' => \$epigenome_name,
   'only_summary!'    => \$only_summary,
);

Bio::EnsEMBL::Registry->load_all($registry);

use Bio::EnsEMBL::Utils::Logger;
my $logger = Bio::EnsEMBL::Utils::Logger->new();

my $ontology_term_adaptor      = Bio::EnsEMBL::Registry->get_adaptor( 'Multi', 'Ontology', 'OntologyTerm' );
my $regulatory_feature_adaptor = Bio::EnsEMBL::Registry->get_adaptor( $species, 'Funcgen', 'RegulatoryFeature' );
my $epigenome_adaptor          = Bio::EnsEMBL::Registry->get_adaptor( $species, 'Funcgen', 'Epigenome' );

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
  $serialiser_callback = create_file_and_serialiser_for_summary($parameters);
} else {
  $serialiser_callback = create_file_and_serialiser($parameters);
}

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

# Avoid error 99 mysql problem. Otherwise the serialisers will connect and 
# disconnect to the ontology (gff serialiser) or core (bed serialiser for 
# fetching ucsc synonyms) database and eventually fail with aforementioned 
# error.
#
# Especially tricky, because the core database connection is never used in the
# script, but heavily used when serialising all features to bed.
#
my $core_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'Core' );
$core_adaptor->dbc->disconnect_when_inactive(0);

# Try to prevent 
#
# DBD::mysql::st execute failed: Lost connection to MySQL server during query 
# at /nfs/users/nfs_m/mn1/work_dir_ftp/lib/ensembl/modules/Bio/EnsEMBL/DBSQL/SliceAdaptor.pm 
# line 825. DBD::mysql::st execute failed: Lost connection to MySQL server 
# during query at /nfs/users/nfs_m/mn1/work_dir_ftp/lib/ensembl/modules/Bio/EnsEMBL/DBSQL/SliceAdaptor.pm 
# line 825.
#
$core_adaptor->dbc->reconnect_when_lost(1);
$ontology_term_adaptor->db->dbc->disconnect_when_inactive(0);
$funcgen_adaptor->dbc->disconnect_when_inactive(0);

while ($exported_something) {

  $exported_something = undef;

  $helper->execute_no_return(
    -SQL      => 'select regulatory_feature_id from regulatory_feature join regulatory_build using (regulatory_build_id) where is_current=1 and regulatory_feature_id > ? order by regulatory_feature_id limit ?',
    -PARAMS => [ $last_id, $batch_size ],
    -CALLBACK => $serialiser_callback,
  );
}

$logger->info("Export done.\n");

sub create_file_and_serialiser_for_summary {

  my $param = shift;

  my $logger                     = $param->{logger};
  my $regulatory_feature_adaptor = $param->{regulatory_feature_adaptor};
  my $exported_something         = $param->{exported_something};
  
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
  return $regulatory_activity_serialiser;
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
  return $regulatory_activity_serialiser;
}
