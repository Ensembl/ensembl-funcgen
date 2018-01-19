#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

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

__END__

>perl metadata_4.pl

AP1 has 0 DBentries
AP1 has 2 associated feature types
        Cfos is associated to ENSG00000170345
        Cjun is associated to ENSG00000177606
