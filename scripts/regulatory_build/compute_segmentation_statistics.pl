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

  compute_segmentation_statistics.pl \
    --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb7.pm \
    --species homo_sapiens

rm -rf /hps/nobackup2/production/ensembl/mnuhn/segmentation_files_while_maintenance_is_ongoing/segmentation_statistics ; compute_segmentation_statistics.pl     --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb7.pm     --species homo_sapiens     --db_file_path /hps/nobackup2/production/ensembl/mnuhn/segmentation_files_while_maintenance_is_ongoing/dbfiles     --tempdir /hps/nobackup2/production/ensembl/mnuhn/segmentation_files_while_maintenance_is_ongoing/segmentation_statistics

=cut

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use Bio::EnsEMBL::Utils::Logger;
use Bio::EnsEMBL::Funcgen::SegmentationStatistic;

my %options;
GetOptions (
    \%options,
    "species|s=s",
    "registry|r=s",
 );

use Hash::Util qw( lock_keys );
lock_keys( %options );

my $species      = $options{'species'};
my $registry     = $options{'registry'};

Bio::EnsEMBL::Registry->load_all($registry);
my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

my $segmentation_adaptor           = $funcgen_adaptor->get_SegmentationAdaptor;
my $segmentation_statistic_adaptor = $funcgen_adaptor->get_SegmentationStatisticAdaptor;

my $all_segmentations = $segmentation_adaptor->fetch_all;

my @segmentation_assignments = qw(
  ctcf
  dead
  distal
  gene
  poised
  proximal
  repressed
  tss
  weak
);

foreach my $assignment (@segmentation_assignments) {
    $logger->info("Computing overall statistics for $assignment");
    
    my $stats = compute_segmentation_statistics_for_assignment(
      $assignment
    );

    $segmentation_statistic_adaptor->_delete_possibly_existing_statistics($stats);
    $segmentation_statistic_adaptor->store($stats);
}

foreach my $assignment (@segmentation_assignments) {

  foreach my $segmentation (@$all_segmentations) {

    $logger->info("Computing statistics for " . $segmentation->name . " and $assignment");
    
    my $stats = compute_segmentation_statistics_for_segmentation_assignment(
      $segmentation,
      $assignment
    );

    $segmentation_statistic_adaptor->_delete_possibly_existing_statistics($stats);
    $segmentation_statistic_adaptor->store($stats);
  }
}

$logger->finish_log;
exit(0);

sub compute_segmentation_statistics_for_assignment {

  my $assignment = shift;
  
  my $fetch_feature_lengths_by_assignment_and_segmentation_sql
    = q( 
      select 
        length
      from 
        segmentation 
        join segmentation_feature using (segmentation_id)
      where
        assignment = ?
    );

  my $sth 
    = $funcgen_adaptor
        ->dbc
        ->prepare(
          $fetch_feature_lengths_by_assignment_and_segmentation_sql
        );

  $sth->execute($assignment);

  my $statistics = _compute_segmentation_statistics_from_sth($sth);
  
  my $stats = [ 
    map { 
      $_->label( $assignment );
      $_ 
    } 
    @{ 
      $statistics
    }
  ];
  return $stats;
}

sub compute_segmentation_statistics_for_segmentation_assignment {

  my $segmentation = shift;
  my $assignment   = shift;
  
  my $fetch_feature_lengths_by_assignment_and_segmentation_sql
    = q( 
      select 
        length
      from 
        segmentation 
        join segmentation_feature using (segmentation_id)
      where
        segmentation.name = ?
        and assignment    = ?
    );

  my $sth 
    = $funcgen_adaptor
        ->dbc
        ->prepare(
          $fetch_feature_lengths_by_assignment_and_segmentation_sql
        );

  $sth
    ->execute(
      $segmentation->name,
      $assignment,
    );

  my $statistics = _compute_segmentation_statistics_from_sth($sth);
  
  my $stats = [ 
    map { 
      $_->label           ( $assignment   );
      $_->set_Segmentation( $segmentation );
      $_ 
    } 
    @{ 
      $statistics
    }
  ];
  return $stats;
}

sub _compute_segmentation_statistics_from_sth {
  my $sth = shift;
  my $lengths = [
    map {
      # Flatten the array returned from the DBI
      #
      @$_ 
    } @{
      $sth->fetchall_arrayref
    }
  ];
  my $statistics = analyse_data( $lengths );
  undef $lengths;
  return $statistics;
}

sub analyse_data {

  my $segmentation_lengths = shift;

  my @segmentation_statistics;
  
  use Statistics::Descriptive;
  my $stat = Statistics::Descriptive::Full->new();
  $stat->add_data(@$segmentation_lengths);
  
  my $q0 = $stat->quantile(0);
  my $q1 = $stat->quantile(1);
  my $q2 = $stat->quantile(2);
  my $q3 = $stat->quantile(3);
  my $q4 = $stat->quantile(4);
  
  if (! defined $q0) { $q0 = 0; }
  if (! defined $q1) { $q1 = 0; }
  if (! defined $q2) { $q2 = 0; }
  if (! defined $q3) { $q3 = 0; }
  if (! defined $q4) { $q4 = 0; }

  $logger->info("\n");
  $logger->info("    Q0 " . $q0 . "\n");
  $logger->info("    Q1 " . $q1 . "\n");
  $logger->info("    Q2 " . $q2 . "\n");
  $logger->info("    Q3 " . $q3 . "\n");
  $logger->info("    Q4 " . $q4 . "\n");

  push 
    @segmentation_statistics, 
    Bio::EnsEMBL::Funcgen::SegmentationStatistic->new(
      -statistic => 'length_q0',
      -value     => $q0,
    ),
    Bio::EnsEMBL::Funcgen::SegmentationStatistic->new(
      -statistic => 'length_q1',
      -value     => $q1,
    ),
    Bio::EnsEMBL::Funcgen::SegmentationStatistic->new(
      -statistic => 'length_q2',
      -value     => $q2,
    ),
    Bio::EnsEMBL::Funcgen::SegmentationStatistic->new(
      -statistic => 'length_q3',
      -value     => $q3,
    ),
    Bio::EnsEMBL::Funcgen::SegmentationStatistic->new(
      -statistic => 'length_q4',
      -value     => $q4,
    ),
  ;

  my $skewness = $stat->skewness;
  my $kurtosis = $stat->kurtosis;
  
  if (! defined $skewness) { $skewness = 'null'; }
  if (! defined $kurtosis) { $kurtosis = 'null'; }
  
  $logger->info("\n");
  $logger->info("    Skewness: " . $skewness . "\n");
  $logger->info("    Kurtosis: " . $kurtosis . "\n");
  
  $logger->info("\n");
  $logger->info("    Sum:   " . $stat->sum   . "\n");
  $logger->info("    Count: " . $stat->count . "\n");
  $logger->info("    Mean:  " . $stat->mean  . "\n");

  push 
    @segmentation_statistics, 
    Bio::EnsEMBL::Funcgen::SegmentationStatistic->new(
      -statistic => 'skewness',
      -value     => $skewness,
    ),
    Bio::EnsEMBL::Funcgen::SegmentationStatistic->new(
      -statistic => 'kurtosis',
      -value     => $kurtosis,
    ),
    Bio::EnsEMBL::Funcgen::SegmentationStatistic->new(
      -statistic => 'sum',
      -value     => $stat->sum,
    ),
    Bio::EnsEMBL::Funcgen::SegmentationStatistic->new(
      -statistic => 'count',
      -value     => $stat->count,
    ),
    Bio::EnsEMBL::Funcgen::SegmentationStatistic->new(
      -statistic => 'mean',
      -value     => $stat->mean,
    ),
  ;

  return \@segmentation_statistics;
}

