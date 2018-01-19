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
use Bio::EnsEMBL::Registry;

Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

my $gene_adaptor          = Bio::EnsEMBL::Registry->get_adaptor("mus musculus", "Core",    "Gene");
my $transcript_adaptor    = Bio::EnsEMBL::Registry->get_adaptor("mus musculus", "Core",    "Transcript");
my $slice_adaptor         = Bio::EnsEMBL::Registry->get_adaptor("mus musculus", "Core",    "Slice");

my $array_adaptor         = Bio::EnsEMBL::Registry->get_adaptor("mus musculus", "Funcgen", "Array");
my $probe_feature_adaptor = Bio::EnsEMBL::Registry->get_adaptor("mus musculus", "Funcgen", "ProbeFeature");
my $probe_set_adaptor     = Bio::EnsEMBL::Registry->get_adaptor("mus musculus", "Funcgen", "ProbeSet");

#####################################################
# Get the slice for gene Cntnap1 with 100000bp flanks
# Get the MOE430A array
#####################################################

# Fetch first gene that goes by Cntnap1
my ($gene) = @{$gene_adaptor->fetch_all_by_external_name('Cntnap1')};
my $slice = $gene->feature_Slice->expand(10000, 10000);
my $array = $array_adaptor->fetch_by_name_vendor('MG-U74B', 'AFFY');


# Get the locations to which this probe has mapped, "probe features" on this 
# slice and array
#
my $probe_features = $probe_feature_adaptor->fetch_all_by_Slice_Array($slice,$array);

print "\n".scalar(@$probe_features)." ProbeFeatures for Array ".$array->name." on slice ".$slice->name."\n\n";

foreach my $current_probe_feature (@$probe_features) {

  my $probeset   = $current_probe_feature->probeset->name;
  my $probename  = $current_probe_feature->probe->get_probename($array->name);
  my $slice_name = $current_probe_feature->feature_Slice->name;

  print "\tProbe $probeset $probename is aligned on $slice_name with cigar string ".$current_probe_feature->cigar_string."\n";
}

############################################
# Get all ProbeFeatures annotated on the probesets annotated to Transcripts
# Note: ProbeFeature level DBEntries are used to build Probe/ProbeSet
# level transcript annotatations.
print "\nProbeFeatures xref'd to Transcripts on slice ".$slice->name."\n";
my @genes = @{ $slice->get_all_Genes() };

foreach my $gene (@genes) {

  print "\n",$gene->external_name," ",$gene->display_id(),"\n";

  foreach my $transcript (@{$gene->get_all_Transcripts}) {

    print "  ", $transcript->external_name, " ", $transcript->display_id(), " has the following probes matches:\n";

    my $probe_feature_on_transcript = $probe_feature_adaptor->fetch_all_by_external_name($transcript->stable_id);
    
    PROBE_FEATURE:
    foreach my $current_probe_feature (@$probe_feature_on_transcript) {

      # Now filter ProbeFeatures which are on our array
      my $on_our_array = 0;
      foreach my $array_having_current_probe_feature (@{$current_probe_feature->probe->get_all_Arrays}) {
        if($array_having_current_probe_feature->name eq $array->name) {
          $on_our_array = 1;
          last;
        }
      }
      next PROBE_FEATURE if ! $on_our_array;

      my $probeset = $current_probe_feature->probeset->name;

      # Get the DBEntry information
      # Filter for transcript by passing optional argument
      foreach my $current_db_entry (@{$current_probe_feature->get_all_Transcript_DBEntries($transcript)}) {

        # Probe can have different names on different arrays
        # So have to specify which array we are dealing with
        #
        my $complete_name = $current_probe_feature->probe->get_complete_name($array->name);
        my $slice_name = $current_probe_feature->feature_Slice->name;
        my $dbe_info = $current_db_entry->linkage_annotation;

        print "\t$complete_name $slice_name\t\tHits $dbe_info\n";
      }
    }
  }
}

# Get All the probesets that have been mapped to a transcript
#
my $demo_transcript_stable_id = 'ENSMUST00000103109';

print "\nAll probesets annotated on $demo_transcript_stable_id:\n";

my $transcript = $transcript_adaptor->fetch_by_stable_id($demo_transcript_stable_id);

if (! defined $transcript) {
  die("Can't find transcript with stable id ${demo_transcript_stable_id}!");
}

my $probesets = $probe_set_adaptor->fetch_all_by_external_name($transcript->stable_id);

foreach my $probeset (@$probesets) {

  my $arrays_string = join(', ', (map $_->name, @{$probeset->get_all_Arrays}));
  my $dbe_info;

  #Now get linkage_annotation
  foreach my $current_db_entry (@{$probeset->get_all_Transcript_DBEntries($transcript)}) {
    $dbe_info = $current_db_entry->linkage_annotation;
  }
  print "\t".$probeset->name." on arrays ".$arrays_string." with Probe hits $dbe_info\n";
}

#############################################################
# Get all probesets annotated to transcripts for a given gene
#############################################################

print "\nAll probesets associated with transcript of Cntnap1:\n";

foreach my $probe_set (@{$probe_set_adaptor->fetch_all_by_linked_transcript_Gene($gene)}) {
  print $probe_set->name."\n";
}

############################################
# Probeset centric access
############################################

print "\nProbeset centric access. All transcripts annotated to  ProbeSet 96567_at";

my $probeset = $probe_set_adaptor->fetch_by_array_probeset_name('MG-U74A', '96567_at');

my $transcript_db_entries = $probeset->get_all_Transcript_DBEntries;

foreach my $current_db_entry (@$transcript_db_entries) {

  # Fetch transcript using primary ID of DBEntry object
  my $transcript = $transcript_adaptor->fetch_by_stable_id($current_db_entry->primary_id);
  my $dbe_info = $current_db_entry->linkage_annotation;

  # Grab the gene using the transcript stable ID stored in the DBEntry primary ID
  my $gene_sid = $gene_adaptor->fetch_by_transcript_stable_id($current_db_entry->primary_id)->stable_id;

  print "\t" 
    . $current_db_entry->db_display_name 
    . ' ' . $current_db_entry->primary_id 
    . ' ' . "($gene_sid) at ".$transcript->feature_Slice->name 
    . ' ' . " with Probe hits $dbe_info\n";
}
