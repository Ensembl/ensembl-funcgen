use strict;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db
  (
   -host => 'ensembldb.ensembl.org',
   -user => 'anonymous',
  );


#Grab the adaptors
my $efg_db        = $registry->get_DBAdaptor('Human', 'funcgen');
my $array_adaptor = $efg_db->get_ArrayAdaptor;
#my $ac_adaptor    = $efg_db->get_ArrayChipAdaptor;

#Grab the Nimblegen array
my $array         = $array_adaptor->fetch_by_name_vendor
  ('2005-05-10_HG17Tiling_Set', 'NIMBLEGEN');

#Grab the ArrayChips
my @array_chips   = @{$array->get_ArrayChips};
#my @array_chips   = @{$ac_adaptor->fetch_all_by_array_id($array->dbID)};

#Print some ArrayChip info
foreach my $ac ( @array_chips ){
  print "ArrayChip:".$ac->name."\tDesignID:".$ac->design_id."\n";
}
