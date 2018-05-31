#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::ProbeMapping;
use Getopt::Long;




my $registry;
my $species;

GetOptions (
   'registry=s'            => \$registry,
   'species=s'             => \$species,
);

use Bio::EnsEMBL::Utils::Logger;
my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

Bio::EnsEMBL::Registry->load_all($registry);

my $meta_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'MetaContainer');

$logger->info("Fetching assembly for $species\n");
my $assembly = $meta_adaptor->single_value_by_key('assembly.default',0);
$logger->info("Done.\n");

$logger->info("Fetching gene build for $species\n");
my $gene_build = $meta_adaptor->single_value_by_key('genebuild.last_geneset_update',0);
$logger->info("Done.\n");

$logger->info("Fetching release version for $species\n");
my $release_version = $meta_adaptor->get_schema_version();
$logger->info("Done.\n");

my $release_date_time = localtime;



$logger->info("Insert meta data in probe_mapping table for $species\n");

my $probe_mapping_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'ProbeMapping');

my $remmove_probemapping = $probe_mapping_adaptor->fetch_single_object;
if ($remmove_probemapping){
	$probe_mapping_adaptor->_delete($remmove_probemapping);
}

my $probe_mapping = Bio::EnsEMBL::Funcgen::ProbeMapping->new(

    -assembly  			 => $assembly,
    -gene_build_version  => $gene_build,
    -release_version     => $release_version,
    -release_date        => $release_date_time,

);

$probe_mapping_adaptor->store($probe_mapping);

$logger->info("Done.\n");
