#!/usr/bin/env perl

use strict;
use Bio::EnsEMBL::Registry;
use File::Spec;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::DBSQL::DBConnection;
use Getopt::Long;

my $registry;
my $species;
my $gff_file;
my $bed_file;
my $id_file;

GetOptions (
   'registry=s' => \$registry,
   'species=s'  => \$species,
   'gff_file=s' => \$gff_file,
   'bed_file=s' => \$bed_file,
   'ids=s'      => \$id_file,
);

Bio::EnsEMBL::Registry->load_all($registry);

use Bio::EnsEMBL::Utils::Logger;
my $logger = Bio::EnsEMBL::Utils::Logger->new();

use Bio::EnsEMBL::Funcgen::Ftp::Utils qw( create_file_handle );

my $ontology_term_adaptor = Bio::EnsEMBL::Registry->get_adaptor( 'Multi', 'Ontology', 'OntologyTerm' );

my $gff_fh = create_file_handle($gff_file);

use Bio::EnsEMBL::Utils::IO::GFFSerializer;
my $gff_serializer = Bio::EnsEMBL::Utils::IO::GFFSerializer->new(
  $ontology_term_adaptor,
  $gff_fh
);

my $bed_fh = create_file_handle($bed_file);

use Bio::EnsEMBL::Utils::IO::BEDSerializer;
my $bed_serializer = Bio::EnsEMBL::Utils::IO::BEDSerializer->new(
  $bed_fh
);

my $annotated_feature_adaptor = Bio::EnsEMBL::Registry->get_adaptor( $species, 'Funcgen', 'AnnotatedFeature' );
my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'Funcgen' );

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
$ontology_term_adaptor->db->dbc->disconnect_when_inactive(0);
$funcgen_adaptor->dbc->disconnect_when_inactive(0);

open my $in, "<", $id_file;

while (my $current_annotated_feature_id = <$in>) {

  chomp($current_annotated_feature_id);
  my $annotated_feature = $annotated_feature_adaptor->fetch_by_dbID($current_annotated_feature_id);
  
  #
  # Check, if the peak extends beyond the seq region. If so, trim to the 
  # length of the sequence region. Generating peaks that extend beyond the
  # seq region is a known bug in SWEmbl and these should be corrected before
  # storing them in the database, but sometimes they still creep through.
  #
  my $feature_end    = $annotated_feature->end;
  my $seq_region_end = $annotated_feature->slice->end;
  
  my $feature_extends_beyond_seq_region_bounds = $feature_end > $seq_region_end;
  
  if ($feature_extends_beyond_seq_region_bounds) {
    $annotated_feature->end($seq_region_end);
  }
  
  # Another known bug: peaks that start at position zero.
  if ($annotated_feature->start == 0) {
    $annotated_feature->start(1);
  }

  eval {
    $gff_serializer->print_feature($annotated_feature);
    $bed_serializer->print_feature($annotated_feature);
  };
  if ($@) {
    use Carp;
    use Data::Dumper;
    
    $Data::Dumper::Sortkeys = 1;
    $Data::Dumper::Maxdepth = 3;
    
    confess(
      "Unable to serialise feature! dbid:"
      . $annotated_feature->dbID
      . "\n"
      . Dumper($annotated_feature)
      . "\n"
      . Dumper($annotated_feature->summary_as_hash)
      . "\n"
      . "The error was:" . $@
    );
  }
  # This prevents memory leaks.
  undef %$annotated_feature;
}

$in->close();

$logger->info("Export done.\n");

$gff_fh->close;
$bed_fh->close;
