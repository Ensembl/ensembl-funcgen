#!/usr/bin/env perl

use strict;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Logger;

my %options;

GetOptions (
  \%options,
  "url=s",
  "registry|r=s",
  "species=s",
  "chromosome=s",
  "stable_id_prefix=s",
  "outfile=s",
);

my $species          = $options{species};
my $stable_id_prefix = $options{stable_id_prefix};
my $outfile          = $options{outfile};
my $registry = $options{'registry'};

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

$logger->info("species          = $species\n");
$logger->info("stable_id_prefix = $stable_id_prefix\n");
$logger->info("outfile          = $outfile\n");
$logger->info("registry         = $registry\n");

Bio::EnsEMBL::Registry->load_all($registry);

my $regulatory_build_adaptor   = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'RegulatoryBuild');
my $regulatory_feature_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'RegulatoryFeature');

use File::Basename qw( dirname );

my $output_directory = dirname($outfile);

if (! -e $output_directory) {
    use File::Path qw( mkpath );
    mkpath( $output_directory );
}

open my $out_fh, '>' . $outfile or die("Can't open file ${outfile}!");

my $regulatory_feature_iterator = $regulatory_feature_adaptor->fetch_Iterator;

while (my $current_regulatory_feature = $regulatory_feature_iterator->next) {

  my $bed_string = join("\t", (
    $current_regulatory_feature->seq_region_name,
    $current_regulatory_feature->bound_start,
    $current_regulatory_feature->bound_end,
    $current_regulatory_feature->feature_type->name,
    remove_stable_id_prefix($stable_id_prefix, $current_regulatory_feature->stable_id),
    $current_regulatory_feature->dbID
  ));
  $out_fh->print($bed_string . "\n");
}
$out_fh->close;

sub remove_stable_id_prefix {

  my $stable_id_prefix = shift;
  my $stable_id = shift;
  
  $stable_id =~ s/^$stable_id_prefix//;
  return $stable_id;
}

$logger->finish_log;