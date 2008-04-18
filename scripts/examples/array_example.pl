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

#Grab all the arrays
my @array         = @{$array_adaptor->fetch_all};

#Print some array info
foreach my $array ( @array ){
  print "\nArray:\t".$array->name."\n";
  print "Type:\t".$array->type."\n";
  print "Vendor:\t".$array->vendor."\n";
}
