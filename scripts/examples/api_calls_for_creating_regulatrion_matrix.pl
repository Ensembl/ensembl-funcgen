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
  
  api_calls_for_creating_regulatrion_matrix.pl

=cut


use strict;
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use Carp;

use strict;

Bio::EnsEMBL::Registry->load_registry_from_db(
  -host => 'ensembldb.ensembl.org',
  -user => 'anonymous',
);

my $funcgen_dba = Bio::EnsEMBL::Registry->get_DBAdaptor('Human', 'Funcgen');

my @feature_type_classes = (
  'Transcription Factor',
  'Open Chromatin',
  'Histone',
  'Polymerase',
);

# The columns of the matrix

my $epigenome_adaptor = $funcgen_dba->get_EpigenomeAdaptor;

foreach my $class (@feature_type_classes) {

  print "Class: $class\n";

  my $epigenomes
    = $epigenome_adaptor->fetch_all_having_PeakCalling_by_class($class);
  
  foreach my $epigenome (@$epigenomes) {
    print "  - " . $epigenome->display_label . "\n";
  }
}

# The rows of the matrix

my $feature_type_adaptor = $funcgen_dba->get_FeatureTypeAdaptor;

foreach my $class (@feature_type_classes) {

  print "Class: $class\n";

  my $feature_types 
    = $feature_type_adaptor->fetch_all_having_PeakCalling_by_class($class);
  
  foreach my $feature_type (@$feature_types) {
    print "  - " . $feature_type->name . "\n";
  }
}

# The entries of the matrix

my $peak_calling_adaptor = $funcgen_dba->get_PeakCallingAdaptor;

foreach my $class (@feature_type_classes) {

  my $epigenomes
    = $epigenome_adaptor->fetch_all_having_PeakCalling_by_class($class);

  my $feature_types 
    = $feature_type_adaptor->fetch_all_having_PeakCalling_by_class($class);

  my $num_epigenomes    = @$epigenomes;
  my $num_feature_types = @$feature_types;
  
  print "\n\nClass: $class has $num_feature_types feature types and $num_epigenomes epigenomes.\n\n\n";
  
  foreach my $epigenome (@$epigenomes) {
  
    my $epigenome_name = $epigenome->display_label;
    my $num_peak_calling_found = 0;
    
    foreach my $feature_type (@$feature_types) {
    
      my $peak_callings = $peak_calling_adaptor->fetch_all_by_Epigenome_FeatureType($epigenome, $feature_type);
      my $peak_calling_found = @$peak_callings > 0;
      
      if ($peak_calling_found) {
        $num_peak_calling_found++;
      }
    }
    
    print "The epigenome $epigenome_name has peak callings for $num_peak_calling_found feature types.\n";
    
    if ($num_peak_calling_found == 0) {
      die("This should never happen!");
    }
  }
}


