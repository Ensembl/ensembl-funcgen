#!/usr/bin/env perl

use strict;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Getopt::Long;

# create_chromosome_length_file_for_bigbed_conversion.pl --registry /nfs/users/nfs_m/mn1/work_dir_ftp/lib/ensembl-funcgen/registry.pm --species homo_sapiens

my $registry;
my $species;

# Deliberately not supporting assembly, because fetching all slices from 
# toplevel is not supported from assemblies other than the current one.
#
# my $assembly;

GetOptions (
   'registry=s' => \$registry,
   'species=s'  => \$species,
#    'assembly=s' => \$assembly,
);

Bio::EnsEMBL::Registry->load_all($registry);

my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'Slice');
my $all_slices = $slice_adaptor->fetch_all('toplevel');

foreach my $current_slice (@$all_slices) {

  my $ucsc_name;
  my $seq_region_name = $current_slice->seq_region_name;
  
  my $ucsc_synonyms = $current_slice->get_all_synonyms('UCSC');
  if(@{$ucsc_synonyms}) {
    $ucsc_name = $ucsc_synonyms->[0]->name();
  } else {
#     $ucsc_name = $current_slice->seq_region_name;
    if($current_slice->is_chromosome()) {
      #MT is a special case; it's chrM
      if($seq_region_name eq 'MT' ) {
        $ucsc_name = 'chrM';
      }
      # If it was a ref region add chr onto it (only check if we have an adaptor)
      elsif($current_slice->is_reference) {
        $ucsc_name = 'chr'.$seq_region_name;
      }
    }
  }
  $ucsc_name = $seq_region_name if ! defined $ucsc_name;
  
  print $ucsc_name;
  print "\t";
  print $current_slice->seq_region_length;
  
  print "\n";
}
