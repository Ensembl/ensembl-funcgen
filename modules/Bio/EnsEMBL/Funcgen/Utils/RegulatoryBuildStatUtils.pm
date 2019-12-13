package Bio::EnsEMBL::Funcgen::Utils::RegulatoryBuildStatUtils;
=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut
use strict;
use warnings;
use base qw( Exporter );
use vars qw( @EXPORT_OK );

use constant {

  CTCF            => 'CTCF Binding Site',
  ENHANCER        => 'Enhancer',
  PROMOTER_FLANK  => 'Promoter Flanking Region',
  PROMOTER        => 'Promoter',
  TF              => 'TF binding site',
  OPEN_CHROMATIN  => 'Open chromatin',

};

use constant 
  REGULATORY_FEATURE_TYPES => (
    CTCF,
    ENHANCER,
    PROMOTER_FLANK,
    PROMOTER,
    TF,
    OPEN_CHROMATIN
  ),
;

@EXPORT_OK = qw(
  range_register_regulatory_features
  REGULATORY_FEATURE_TYPES
  CTCF
  ENHANCER
  PROMOTER_FLANK
  PROMOTER
  TF
  OPEN_CHROMATIN
);

sub range_register_regulatory_features {

  my $param = shift;
  my $species = $param->{species};
  my $max     = $param->{max};
  
  my $stop_after_maximum_reached = undef;
  
  if (defined $max) {
    $stop_after_maximum_reached = 1;
  }
  
  use Bio::EnsEMBL::Mapper::RangeRegistry;
  my $range_registry = Bio::EnsEMBL::Mapper::RangeRegistry->new();
  
  my $core_adaptor    = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'core');
  my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
  
  my $regulatory_build_adaptor   = $funcgen_adaptor->get_RegulatoryBuildAdaptor;
  my $regulatory_feature_adaptor = $funcgen_adaptor->get_RegulatoryFeatureAdaptor;
  my $regulatory_build = $regulatory_build_adaptor->fetch_current_regulatory_build;

  my $iterator = $regulatory_feature_adaptor->fetch_Iterator_by_RegulatoryBuild($regulatory_build);
  
  my %by_feature_type;
  my %total_overlap_size_by_feature_type;

  foreach my $feature_type (REGULATORY_FEATURE_TYPES) {
    $by_feature_type{$feature_type} = Bio::EnsEMBL::Mapper::RangeRegistry->new();
    $total_overlap_size_by_feature_type{$feature_type} = 0;
  }

  use Hash::Util qw( lock_keys );
  lock_keys(%by_feature_type);

  my $i = 0;
  
  REGULATORY_FEATURE:
  while ($iterator->has_next) {

    my $current_regulatory_feature = $iterator->next;
    
    my $feature_type = $current_regulatory_feature->feature_type->name;
    
    $by_feature_type{$feature_type}->check_and_register(
        $current_regulatory_feature->slice->seq_region_name,
#         $current_regulatory_feature->start,
        $current_regulatory_feature->bound_seq_region_start,
#         $current_regulatory_feature->end,
        $current_regulatory_feature->bound_seq_region_end,
#         $current_regulatory_feature->start,
        $current_regulatory_feature->bound_seq_region_start,
#         $current_regulatory_feature->end,
        $current_regulatory_feature->bound_seq_region_end,
    );
    
    $range_registry->check_and_register(
        $current_regulatory_feature->slice->seq_region_name,
#         $current_regulatory_feature->start,
        $current_regulatory_feature->bound_seq_region_start,
#         $current_regulatory_feature->end,
        $current_regulatory_feature->bound_seq_region_end,
#         $current_regulatory_feature->start,
        $current_regulatory_feature->bound_seq_region_start,
#         $current_regulatory_feature->end,
        $current_regulatory_feature->bound_seq_region_end,
    );
    $i++;
    if ($stop_after_maximum_reached && $i == $max) {
      last REGULATORY_FEATURE
    }
  }
  return ( $range_registry, \%by_feature_type );
}

1;
