#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db
  (
   -host => 'ensembldb.ensembl.org',
   -user => 'anonymous',
   -db_version => 67,
   -verbose => 1,
  );


#Grab the adaptors
my $efg_db        = $registry->get_DBAdaptor('homo_sapiens', 'funcgen');#'Human', 'funcgen');
my $array_adaptor = $efg_db->get_ArrayAdaptor;

#Grab all the arrays
my @array         = @{$array_adaptor->fetch_all};

#Print some array info
foreach my $array ( @array ){
  print "\nArray:\t".$array->name."\n";
  print "Type:\t".$array->type."\n";
  print "Vendor:\t".$array->vendor."\n";
}
