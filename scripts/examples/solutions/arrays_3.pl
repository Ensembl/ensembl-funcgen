#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
