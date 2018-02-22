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

  perl scripts/regulatory_build/compute_genome_coverage.pl --registry /homes/mnuhn/work_dir_regbuild_script/lib/ensembl-funcgen/registry.pm --species homo_sapiens --tempdir foobar
  perl scripts/regulatory_build/compute_genome_coverage.pl --registry /homes/mnuhn/work_dir_regbuild_script/lib/ensembl-funcgen/registry.pm --species mus_musculus --tempdir foobar

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

use Bio::EnsEMBL::Funcgen::Utils::RegulatoryBuildStatUtils qw ( 
  range_register_regulatory_features 
  REGULATORY_FEATURE_TYPES
);

(
  my $range_registry,
  my $by_feature_type,
)
  = range_register_regulatory_features({
    species => $species
  });

my $genome_size = 0;
my $total_overlap_size = 0;

my %total_overlap_size_by_feature_type;

my $slice_adaptor = $core_adaptor->get_SliceAdaptor;
my $karyotype_slice = $slice_adaptor->fetch_all_karyotype;

foreach my $current_karyotype_slice (@$karyotype_slice) {

  my $seq_region_name = $current_karyotype_slice->seq_region_name;
  my $length          = $current_karyotype_slice->length;
  
  $genome_size += $length;
  
  my $overlap_size = $range_registry->overlap_size($seq_region_name, 0, 1 + $length);
  $total_overlap_size += $overlap_size;
  
  foreach my $feature_type (REGULATORY_FEATURE_TYPES) {
    $total_overlap_size_by_feature_type{$feature_type} += $by_feature_type->{$feature_type}->overlap_size($seq_region_name, 0, 1 + $length);
  }
}

my %statistics;

my $percent_overlap = 100 * $total_overlap_size / $genome_size;

$statistics{genome_size}        = $genome_size;
$statistics{total_overlap_size} = $total_overlap_size;
$statistics{percent_overlap}    = $percent_overlap;

foreach my $feature_type (REGULATORY_FEATURE_TYPES) {

  my $total_overlap_size = $total_overlap_size_by_feature_type{$feature_type};
  $percent_overlap = 100 * $total_overlap_size / $genome_size;
  
  my %feature_type_statistics;
  
  $feature_type_statistics{total_overlap_size} = $total_overlap_size;
  $feature_type_statistics{percent_overlap}    = $percent_overlap;
  
  $statistics{$feature_type} = \%feature_type_statistics;
}

print "Species:,$species\n";
print "Genome size:,$statistics{genome_size}\n";
print "\n";
print "Feature type:,Total bp covered,Percent of genome covered\n";
print "All,$statistics{total_overlap_size},$statistics{percent_overlap}\n";

foreach my $feature_type (REGULATORY_FEATURE_TYPES) {

  print "$feature_type,$statistics{$feature_type}{total_overlap_size},$statistics{$feature_type}{percent_overlap}\n";

}

$logger->finish_log;
