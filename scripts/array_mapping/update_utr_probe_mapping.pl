#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Getopt::Long;
use Hash::Util qw( lock_hash );

my $registry;
my $species;
my $unannotated_utrs_file;

GetOptions (
   'registry=s'                       => \$registry,
   'species=s'                        => \$species,
   'unannotated_utrs=s'               => \$unannotated_utrs_file,
);

use Bio::EnsEMBL::Utils::Logger;
my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

Bio::EnsEMBL::Registry->load_all($registry);

$logger->info("Fetching utr data from $unannotated_utrs_file\n");
use Bio::EnsEMBL::Funcgen::Parsers::DataDumper;
my $unannotated_utrs = Bio::EnsEMBL::Funcgen::Parsers::DataDumper->new->load_first_item_from_data_dump_file($unannotated_utrs_file);
lock_hash(%$unannotated_utrs);
if (! defined $unannotated_utrs){
	die("unannotated_utrs not defined");
}
$logger->info("Done\n");

$logger->info("Updating probe_mapping table\n");
my $probe_mapping_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'ProbeMapping');
my $probe_mapping = $probe_mapping_adaptor->fetch_single_object;

$probe_mapping->three_prime_utr($unannotated_utrs->{3});
$probe_mapping->five_prime_utr($unannotated_utrs->{5});

$probe_mapping_adaptor->update($probe_mapping);


$logger->info("Done\n");

