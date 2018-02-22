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
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=head1 NAME

  A script to compute the number of bases covered by regulatory features both in total and by feature type

  perl scripts/regulatory_build/compute_enhancer_coverage.pl --registry /homes/mnuhn/work_dir_regbuild_script/lib/ensembl-funcgen/registry.pm --species homo_sapiens --tempdir foobar
  perl scripts/regulatory_build/compute_enhancer_coverage.pl --registry /homes/mnuhn/work_dir_regbuild_script/lib/ensembl-funcgen/registry.pm --species mus_musculus --tempdir foobar

=cut

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Mapper::RangeRegistry;
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::Slice;

use Bio::EnsEMBL::Mapper::RangeRegistry;

my $range_registry = Bio::EnsEMBL::Mapper::RangeRegistry->new();

use Bio::EnsEMBL::Utils::Logger;

my %options;
GetOptions (
    \%options,
    "species|s=s",
    "registry|r=s",
    "tempdir|t=s",
 );

use Hash::Util qw( lock_keys );
lock_keys( %options );

my $species  = $options{'species'};
my $registry = $options{'registry'};
my $tempdir  = $options{'tempdir'};

Bio::EnsEMBL::Registry->load_all($registry);

my $core_adaptor    = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'core');
my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

my $feature_set_adaptor = $funcgen_adaptor->get_FeatureSetAdaptor;
my $external_feature_adaptor = $funcgen_adaptor->get_ExternalFeatureAdaptor;

my @experimentally_verified_enhancer_feature_set_names = (
  'FANTOM',
  'VISTA enhancer set',
);

my @experimentally_verified_enhancer_feature_sets
  = grep { 
    defined $_
  } map {
    $feature_set_adaptor->fetch_by_name($_);
  } @experimentally_verified_enhancer_feature_set_names;

my $slice_adaptor = $core_adaptor->get_SliceAdaptor;
my $karyotype_slice = $slice_adaptor->fetch_all_karyotype;

use Bio::EnsEMBL::Funcgen::Utils::RegulatoryBuildStatUtils qw ( 
  range_register_regulatory_features 
  REGULATORY_FEATURE_TYPES
  ENHANCER
);

$logger->info("Range registering regulatory features\n");

(
  my $range_registry,
  my $by_feature_type,
)
  = range_register_regulatory_features({
    species => $species,
#     max     => 1000,
  });

my $registered_enhancers = $by_feature_type->{ENHANCER()};

foreach my $experimentally_verified_enhancer_feature_set (@experimentally_verified_enhancer_feature_sets) {

  my $current_feature_set = $experimentally_verified_enhancer_feature_set->name;
  
  my $total_enhancers_checked       = 0;
  my $num_enhancers_overlapping     = 0;
  my $num_enhancers_not_overlapping = 0;

  $logger->info("Computing coverage of $current_feature_set\n");

  $logger->info("Iterating over karyotype slices\n");

  foreach my $current_karyotype_slice (@$karyotype_slice) {

    my $seq_region_name = $current_karyotype_slice->seq_region_name;
    my $length          = $current_karyotype_slice->length;
    
    $logger->info("    Checking chromosome $seq_region_name\n");
    
    my $feature_iterator = $external_feature_adaptor
      ->fetch_Iterator_by_Slice_FeatureSets(
        $current_karyotype_slice, 
        [ $experimentally_verified_enhancer_feature_set ]
      );
    
    FEATURE:
    while ($feature_iterator->has_next) {
    
      # This may or may not actually be an enhancer
      my $current_enhancer_feature = $feature_iterator->next;
      
      # Check it really is an enhancer, sadly enhancers are hodgepodged together
      # with other stuff in the same feature set.
      #
      if ($current_enhancer_feature->feature_type->class ne ENHANCER) {
        next FEATURE;
      }
      
      my $overlap_size = $registered_enhancers->overlap_size(
        $seq_region_name, 
        $current_enhancer_feature->start, 
        $current_enhancer_feature->end
      );
      
      my $overlaps_enhancer_from_regulatory_build = $overlap_size > 0;
      
      $total_enhancers_checked++;
      
      if ($overlaps_enhancer_from_regulatory_build) {
        $num_enhancers_overlapping++;
      } else {
        $num_enhancers_not_overlapping++;
      }
    }
  }

  $logger->info("Done checking\n");

  print <<REPORT

  Summary for: $current_feature_set
  ============
  
  total_enhancers_checked       = $total_enhancers_checked
  num_enhancers_overlapping     = $num_enhancers_overlapping
  num_enhancers_not_overlapping = $num_enhancers_not_overlapping

REPORT
  ;

}

$logger->finish_log;
