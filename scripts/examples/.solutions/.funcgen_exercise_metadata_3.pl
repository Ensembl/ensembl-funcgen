#!/usr/bin/perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

# 3. Feature Type Gene Associations
# A Feature Type corresponds frequently to an antibody used in a ChIP experiment.
# We associate such Feature Types to Ensembl Genes via External References.
# Print the stable id of the Human Ensembl Gene associated to the Feature Type 'CTCF'.
# Hint: use $efg_dba->get_DBEntryAdaptor() and the class Bio::EnsEMBL::Funcgen::DBSQL::DBEntryAdaptor
# Print the name of the Human Feature Type associated to the Ensembl Gene 'ENSG00000125952'

#Grab the eFG adaptor
my $fta = $registry->get_adaptor('Human', 'funcgen', 'featuretype');
my $dbea = $registry->get_adaptor('Human', 'funcgen', 'dbentry');

my $ctcf = $fta->fetch_by_name('CTCF');
my @db_entries = @{$dbea->fetch_all_by_FeatureType($ctcf)};
#We know there is only one (an Ensembl Gene), but there could be more, so we better try and check them all:
foreach my $db_entry (@db_entries){
	print $db_entry->display_id."\t".$db_entry->primary_id."\n";	
}

my @ft_ids = $dbea->list_feature_type_ids_by_extid("ENSG00000125952");
my $ft = $fta->fetch_by_dbID($ft_ids[0]);
print $ft->name."\n";
