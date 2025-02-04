#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2025] EMBL-European Bioinformatics Institute

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

  select_peak_callings_to_use_in_regulatory_build.pl \
    --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb16.pm \
    --species homo_sapiens \
    --dry_run 0

  classify_epigenome_to_segmentation_run.pl \
      --species  homo_sapiens \
      --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb16.pm \
      --partition_by_experimental_group 1

=cut

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use Bio::EnsEMBL::Utils::Logger;

my %options;
GetOptions (
    \%options,
    "species|s=s",
    "registry|r=s",
    "dry_run|d=s",
 );

use Hash::Util qw( lock_keys );
lock_keys( %options );

my $species      = $options{'species'};
my $registry     = $options{'registry'};
my $dry_run      = $options{'dry_run'};

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

$logger->info("registry   = " . $registry   . "\n");
$logger->info("species    = " . $species    . "\n");
$logger->info("dry_run    = " . $dry_run    . "\n");

Bio::EnsEMBL::Registry->load_all($registry);
my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');

my $peak_calling_adaptor = $funcgen_adaptor->get_PeakCallingAdaptor;

my @ids;

#@ids = (20, 21, 22, 23, 24, 25, 26, 27);
@ids = map { $_->dbID } @{$peak_calling_adaptor->fetch_all};

my $num_checked       = 0;
my $num_passed_checks = 0;
my $num_failed_checks = 0;

foreach my $id (@ids) {

  my $peak_calling = $peak_calling_adaptor->fetch_by_dbID($id);

  my $peak_calling_passed = check_peak_calling($peak_calling);
  $num_checked++;

  $peak_calling->used_for_regulatory_build($peak_calling_passed);
  
  if (!$dry_run) {
    $peak_calling_adaptor->update($peak_calling);
  }
  
  if ($peak_calling_passed) {
    $num_passed_checks++;
  } else {
    $num_failed_checks++;
  }
}

$logger->info("Number of peak callings checked: $num_checked\n");
$logger->info("passed: $num_passed_checks\n");
$logger->info("failed: $num_failed_checks\n");

$logger->info("Done\n");
$logger->finish_log;
exit(0);


sub check_peak_calling {

  my $peak_calling = shift;

  $logger->info("\nPeak calling: " . $peak_calling->name . "\n\n");
  
  my $frip_passed         = 0;
  my $chance_passed       = 0;
  my $phantom_peak_passed = 0;
  
  my $broad_peak_calling = $peak_calling->get_Experiment->get_FeatureType->creates_broad_peaks;
  
  my $frip = $peak_calling->get_Frip;
  
  my $frip_threshold = 0.01;
  
  $frip_passed = $frip->frip >= $frip_threshold;
  
  my $chance = $peak_calling->get_Chance;

  $logger->info("  Chance:\n");

  $chance_passed = 1;
  
#   if (!$chance) {
#       $logger->info("    There is no chance entry.\n");
#       $chance_passed = 0
#   }
# 
#   if ($chance) {
#       my $chance_run_failed  = $chance->run_failed;
#       
#       $logger->info("    run_failed  = " . $chance_run_failed  . "\n");
#       
#       # No useful thresholds are given here:
#       # https://github.com/songlab/chance/wiki/CHANCE-Manual
#       
#       my $percent_genome_enriched                      = $chance->percent_genome_enriched;
#       my $control_enrichment_stronger_than_chip_at_bin = $chance->control_enrichment_stronger_than_chip_at_bin;
#       my $first_nonzero_bin_at                         = $chance->first_nonzero_bin_at;
#       
#       $logger->info("    percent_genome_enriched = " . $percent_genome_enriched . "\n");
#       $logger->info("    control_enrichment_stronger_than_chip_at_bin = " . $control_enrichment_stronger_than_chip_at_bin . "\n");
#       $logger->info("    first_nonzero_bin_at = " . $first_nonzero_bin_at . "\n");
#       $logger->info("    \n");
#       
#       $chance_passed 
#         = 
#                !$chance_run_failed
#           && ( $percent_genome_enriched > 10  )
#           && ( $first_nonzero_bin_at    < 10  )
#       ;
#   }

  my $signal_alignment = $peak_calling->get_signal_Alignment;
  my $phantom_peak     = $signal_alignment->get_PhantomPeak;

  $logger->info("  Phantom peaks\n");

  if (!$phantom_peak) {
    $logger->info("    There is no phantom peak entry.\n");
    $phantom_peak_passed = 0;
  }

  if ($phantom_peak) {

      my $nsc                     = $phantom_peak->nsc;
      my $rsc                     = $phantom_peak->rsc;
      my $quality_tag             = $phantom_peak->quality_tag;
      my $phantom_peak_run_failed = $phantom_peak->run_failed;

      $logger->info("    run_failed  = " . $phantom_peak_run_failed  . "\n");
      $logger->info("    nsc         = " . $nsc                      . "\n");
      $logger->info("    rsc         = " . $rsc                      . "\n");
      $logger->info("    quality_tag = " . $quality_tag              . "\n");
      $logger->info("    \n");

      my $broad_peak_nsc_threshold = 1.02;
      my $broad_peak_rsc_threshold = 0.4;

      my $narrow_peak_nsc_threshold = 1.05;
      my $narrow_peak_rsc_threshold = 1;

      if ($broad_peak_calling) {
      
        $phantom_peak_passed 
          = 
                !$phantom_peak_run_failed
            && ( $nsc         > $broad_peak_nsc_threshold )
            && ( $rsc         > $broad_peak_rsc_threshold )
            #&& ( $quality_tag > 0    )
            
      } else {
      
        $phantom_peak_passed 
          = 
                !$phantom_peak_run_failed
            && ( $nsc         > $narrow_peak_nsc_threshold )
            && ( $rsc         > $narrow_peak_rsc_threshold )
            #&& ( $quality_tag > 0    )
      
      }
      
  }

  my $peak_calling_passed = $frip_passed && $chance_passed && $phantom_peak_passed;

  if (!$peak_calling_passed){
	 $peak_calling_passed = 0;
  }

  $logger->info("    Frip passed         = " . $frip_passed         . "\n");
  $logger->info("    Chance passed       = " . $chance_passed       . "\n");
  $logger->info("    Phantom peak passed = " . $phantom_peak_passed . "\n");
  $logger->info("    -----------------------\n");
  $logger->info("    Peak calling passed = " . $peak_calling_passed . "\n");

  return $peak_calling_passed;
}



