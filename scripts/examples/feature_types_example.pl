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
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

my $ftype_adaptor = $registry->get_adaptor('Human', 'funcgen', 'featuretype');
my $fset_adaptor = $registry->get_adaptor('Human', 'funcgen', 'featureset');

#Print all feature sets for Transcription Factors (note that this does not include CTCF, an insulator)
my @tfs = @{$ftype_adaptor->fetch_all_by_class('Transcription Factor')}; 
foreach my $ft (@tfs){
	print "Feature Type: ".$ft->name."\n";
	my @fsets = @{$fset_adaptor->fetch_all_by_FeatureType($ft)};
	print "\t".scalar(@fsets)." Feature Sets available:\n";
	foreach my $fset (@fsets){ 
		print "\t\t".$fset->name."\n"; 
	}
}
