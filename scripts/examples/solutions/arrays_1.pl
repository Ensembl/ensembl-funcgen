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

### Array example code

use strict;
use Bio::EnsEMBL::Registry;

my $reg = 'Bio::EnsEMBL::Registry';

$reg->load_registry_from_db(
							-host => 'ensembldb.ensembl.org',
							-user => 'anonymous',
						   );

my $efg_db = $reg->get_DBAdaptor('Human', 'funcgen');
#my $array1_adaptor = $reg->get_adaptor('Human', 'funcgen', 'Array');
my $array_adaptor = $efg_db->get_ArrayAdaptor;
#These are the same adaptor
#print "$array_adaptor and $array1_adaptor\n";

#List all arrays details
print "Listing all available arrays\n";
my @arrays = @{$array_adaptor->fetch_all};

foreach my $array(@arrays){
  print "\nArray:\t".$array->name."\n";
  print "Type:\t".$array->type."\n";
  print "Vendor:\t".$array->vendor."\n";
}


### ArrayChip Example code
#List all ArrayChip details for the HG17 Nimblegen tiling set

print "\n\nGetting NIMBLEGEN 2005-05-10_HG17Tiling_Set\n";

my $array = $array_adaptor->fetch_by_name_vendor('2005-05-10_HG17Tiling_Set', 'NIMBLEGEN');
my $ac_adaptor = $efg_db->get_ArrayChipAdaptor;

my @achips = @{$ac_adaptor->fetch_all_by_array_id($array->dbID())};
#my @achips = $array->get_ArrayChips;

foreach my $ac(@achips){
  print "ArrayChip:\t".$ac->name."\tDesign ID:\t".$ac->design_id."\n";
}


### Probe example code

my $probe_adaptor = $efg_db->get_ProbeAdaptor;
my $pfeature_adaptor = $efg_db->get_ProbeFeatureAdaptor;

my $probe = $probe_adaptor->fetch_by_array_probe_probeset_name('2005-05-10_HG17Tiling_Set', 'chr22P38797630');

print "Got ".$probe->class." probe ".$probe->get_probename."\n";

my @pfeatures = @{$pfeature_adaptor->fetch_all_by_Probe($probe)};
#could use convinience method
#my @pfeatures = $probe->get_all_ProbeFeatures;

#This should only have one feature as this is a uniquely mapped tiling probe

print "Found ".scalar(@pfeatures)." ProbeFeatures\n";

foreach my $pfeature(@pfeatures){
  print "ProbeFeature found at ".$pfeature->seq_region_name.' '.$pfeature->seq_region_start.' '.$pfeature->seq_region_end."\n";
  print "ProbeFeature Slice name is ".$pfeature->feature_Slice->name."\n";
}


print "End of example code\n";


#1:
#Fetch the Sanger ENCODE3.1.1 array. Print out it's name, type, description and vendor. Hint: use the ArrayAdaptor->fetch by_name method.



my $encode_array = $array_adaptor->fetch_by_name_vendor('ENCODE3.1.1', 'SANGER');
print "encode "  .$encode_array->name.' '.$encode_array->type.' '.$encode_array->vendor.': '.$encode_array->description."\n";
