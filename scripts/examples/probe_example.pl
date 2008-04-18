use strict;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db
  (
   -host => 'ensembldb.ensembl.org',
   -user => 'anonymous',
  );


#Grab the adaptors
my $efg_db           = $registry->get_DBAdaptor('Human', 'funcgen');
my $probe_adaptor    = $efg_db->get_ProbeAdaptor;
my $pfeature_adaptor = $efg_db->get_ProbeFeatureAdaptor;

#Grab a probe from the HG17 array
my $probe            = $probe_adaptor->fetch_by_array_probe_probeset_name
  ('2005-05-10_HG17Tiling_Set', 'chr22P38797630');
print "Got ".$probe->class." probe ".$probe->get_probename."\n";

#Grab the feature associated with this probe
my @pfeatures        = @{$pfeature_adaptor->fetch_all_by_Probe($probe)};
print "\nFound ".scalar(@pfeatures)." ProbeFeatures\n";

#Print some info about the features
foreach my $pfeature ( @pfeatures ){
  print "\nProbeFeature found at:\t".$pfeature->feature_Slice->name."\n";
}
