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
my $min_id;
my $max_id;
my $feature_set_id;

GetOptions (
   'registry=s' => \$registry,
   'species=s'  => \$species,
   'gff_file=s' => \$gff_file,
   'bed_file=s' => \$bed_file,
   'min_id=s'   => \$min_id,
   'max_id=s'   => \$max_id,
   'feature_set_id=s' => \$feature_set_id,
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

my @batch_constraint_list;
if ($min_id) {
  push @batch_constraint_list, qq(annotated_feature_id>=$min_id);
}
if ($max_id) {
  push @batch_constraint_list, qq(annotated_feature_id<=$max_id);
}
if ($feature_set_id) {
  push @batch_constraint_list, qq(feature_set_id=$feature_set_id);
}
my $batch_constraint = join ' and ', @batch_constraint_list;

my $annotated_feature_adaptor = Bio::EnsEMBL::Registry->get_adaptor( $species, 'Funcgen', 'AnnotatedFeature' );
my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'Funcgen' );

my $helper = Bio::EnsEMBL::Utils::SqlHelper->new(
  -DB_CONNECTION => $funcgen_adaptor->dbc
);
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
	$gff_serializer->print_feature($annotated_feature);
	$bed_serializer->print_feature($annotated_feature);
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

$gff_fh->close;
$bed_fh->close;

# Make sure that empty files exist so subsequent analyses don't fail.
# The gff serialiser might not be creating a file, 
#

if (-z $gff_file) {
  `touch $gff_file`;
}
if (-z $bed_file) {
  `touch $bed_file`;
}


