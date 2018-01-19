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
print $ft->name."\tENSG00000125952\n";

__END__


>perl metadata_3.pl

CTCF    ENSG00000102974
Max     ENSG00000125952

