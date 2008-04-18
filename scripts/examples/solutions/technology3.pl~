use strict;
use Bio::EnsEMBL::Registry;

my $reg = 'Bio::EnsEMBL::Registry';

$reg->load_registry_from_db(
							-host => 'ensembldb.ensembl.org',
							-user => 'anonymous',
						   );

my $efg_db = $reg->get_DBAdaptor('Human', 'funcgen');

#3: How many probes are on the ENCODE array?  

my $array_adaptor = $efg_db->get_ArrayAdaptor;
my $encode_array = $array_adaptor->fetch_by_name_vendor('ENCODE3.1.1', 'SANGER');

my $probe_adaptor = $efg_db->get_ProbeAdaptor;

my @eprobes = @{$probe_adaptor->fetch_all_by_Array($encode_array)};

print scalar(@eprobes)." probes on this array\n";

my $cnt = 0;
foreach my $eprobe(@eprobes){
  print $eprobe->get_probename.' is '.$eprobe->length." bases long\n";
  $cnt++;
  last if $cnt == 5;
}
