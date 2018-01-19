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
use warnings;
use Bio::EnsEMBL::Registry;

Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

my $feature_type_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Human', 'Funcgen', 'FeatureType');
my $feature_set_adaptor  = Bio::EnsEMBL::Registry->get_adaptor('Human', 'funcgen', 'FeatureSet');

my @transcription_factor_feature_types =  @{$feature_type_adaptor->fetch_all_by_class('Transcription Factor')}; 

foreach my $current_transcription_factor_feature_type (@transcription_factor_feature_types) {

    print $current_transcription_factor_feature_type->name . "\n";
    
    my $feature_sets = $feature_set_adaptor->fetch_all_by_FeatureType($current_transcription_factor_feature_type);
    
    foreach my $feature_set (@$feature_sets) {
      print "\thas feature set " . $feature_set->name."\n"; 
    }
    print "\n";
}
