### Array example code

use strict;
use Bio::EnsEMBL::Registry;

my $reg = 'Bio::EnsEMBL::Registry';

$reg->load_registry_from_db(
							-host => 'ensembldb.ensembl.org',
							-user => 'anonymous',
						   );

my $efg_db = $reg->get_DBAdaptor('Human', 'funcgen');

#2:
#How many ArrayChips are part of this Array? Get the design IDs for these chips and print them out. Try using the convinience methods in Array.

my $array_adaptor = $efg_db->get_ArrayAdaptor;
my $encode_array = $array_adaptor->fetch_by_name_vendor('ENCODE3.1.1', 'SANGER');


my @achips = @{$encode_array->get_ArrayChips};

my @design_ids = @{$encode_array->get_design_ids};

print 'Found '.scalar(@achips).' ArrayChips with design IDs: '.join("\t", @design_ids)."\n";

