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
use Bio::EnsEMBL::Registry;

Bio::EnsEMBL::Registry->load_registry_from_db(
  -host => 'ensembldb.ensembl.org',
  -user => 'anonymous',
);

#Grab the adaptors
my $efg_db           = Bio::EnsEMBL::Registry->get_DBAdaptor('Human', 'funcgen');
my $probe_adaptor    = $efg_db->get_ProbeAdaptor;
my $pfeature_adaptor = $efg_db->get_ProbeFeatureAdaptor;

#Grab a probe from the HG17 array
my $probe            = $probe_adaptor->fetch_all_by_name('273:421');

print "Got ".$probe->class." probe ".$probe->get_probename."\n";

#Grab the feature associated with this probe
my @pfeatures        = @{$pfeature_adaptor->fetch_all_by_Probe($probe)};
print "\nFound ".scalar(@pfeatures)." ProbeFeatures\n";

#Print some info about the features
foreach my $pfeature ( @pfeatures ){
  print "\nProbeFeature found at:\t".$pfeature->feature_Slice->name."\n";
}
