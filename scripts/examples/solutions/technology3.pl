use strict;
use Bio::EnsEMBL::Registry;

my $reg = 'Bio::EnsEMBL::Registry';

$reg->load_registry_from_db(
							-host => 'ensembldb.ensembl.org',
							-user => 'anonymous',
						   );

my $efg_db = $reg->get_DBAdaptor('Human', 'funcgen');

print "\nQuestion 3: How many probes are on the ENCODE array?\n\n";

#Grab the encode array from the array adaptor
my $array_adaptor = $efg_db->get_ArrayAdaptor;
my $encode_array = $array_adaptor->fetch_by_name_vendor('ENCODE3.1.1', 'SANGER');

#Grab all the probe from the probe adaptor using the encode array
print "Getting probes from the SANGER ENCODE3.1.1 array\n\n";
my $probe_adaptor = $efg_db->get_ProbeAdaptor;
my @eprobes = @{$probe_adaptor->fetch_all_by_Array($encode_array)};

print scalar(@eprobes)." probes on this array\n\n";

print "Here are some details for the first 5:\n\n";

my $cnt = 0;
foreach my $eprobe(@eprobes){
  print $eprobe->get_probename.' is '.$eprobe->length." bases long\n";
  $cnt++;
  last if $cnt == 5;
}

print "\nWow, those PCR probes are a bit longer than your run of the mill oligo probes\n";
