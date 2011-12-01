#!/usr/bin/perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

# 4. Complex Feature Types
# Some feature types correspond to an association between multiple feature types
# See if there are any Human Ensembl Genes associated to the feature type 'AP1'
# List the feature types associated to this feature type. List the ENSEMBL genes associated to these.

#Grab the eFG adaptor
my $fta = $registry->get_adaptor('Human', 'funcgen', 'featuretype');
my $dbea = $registry->get_adaptor('Human', 'funcgen', 'dbentry');

my $ap1 = $fta->fetch_by_name('AP1');

my @db_entries = @{$dbea->fetch_all_by_FeatureType($ap1)};
print $ap1->name." has ".scalar(@db_entries)." DBentries\n";

my @associated_feature_types = @{$fta->fetch_all_by_association($ap1)};
print $ap1->name." has ".scalar(@associated_feature_types)." associated feature types\n";
foreach my $associated_ft (@associated_feature_types) {
	print "\t".$associated_ft->name;
	@db_entries = @{$dbea->fetch_all_by_FeatureType($associated_ft)};
	if(scalar(@db_entries)>0){
		print " is associated to ".$db_entries[0]->primary_id;
	} else {
		print " has no Associated DBEntry";
	}
	print "\n";
}

