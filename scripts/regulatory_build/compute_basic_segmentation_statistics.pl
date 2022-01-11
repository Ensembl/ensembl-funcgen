#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2022] EMBL-European Bioinformatics Institute

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

  This script requires that load_segmentation_files_to_db.pl was been run first.

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
 );

use Hash::Util qw( lock_keys );
lock_keys( %options );

my $species  = $options{'species'};
my $registry = $options{'registry'};

Bio::EnsEMBL::Registry->load_all($registry);
my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

$logger->info("Creating indexes\n");

eval {
  $funcgen_adaptor->dbc->do(
    q~
      create index foo on segmentation_feature(segmentation_id);
    ~
  );
};
eval {
  $funcgen_adaptor->dbc->do(
    q~
      create index foo2 on segmentation_feature(segmentation_id, assignment);
    ~
  );
};

$logger->info("Computing statistics\n");

$funcgen_adaptor->dbc->do(
  q~
  insert into segmentation_statistic 
    (
      segmentation_id, 
      statistic, 
      value
    )
    select 
      segmentation.segmentation_id, 
      'num_epigenomes', 
      count(distinct epigenome_id) 
    from 
      segmentation 
      join segmentation_cell_tables using (segmentation_id) 
    group by 
      segmentation.segmentation_id
  ~
);

$funcgen_adaptor->dbc->do(
  q~
  insert into segmentation_statistic 
    (
      segmentation_id,
      label,
      statistic, 
      value
    )
    select 
      segmentation.segmentation_id, 
      assignment,
      'average_length',
      avg(length)
    from 
      segmentation 
      join segmentation_feature using (segmentation_id) 
    group by 
      segmentation_id, 
      assignment
  ;
  ~
);

$funcgen_adaptor->dbc->do(
  q~
  insert into segmentation_statistic 
    (
      segmentation_id,
      label,
      statistic, 
      value
    )
    select 
      segmentation.segmentation_id, 
      assignment,
      'num_states',
      count(distinct state)
    from 
      segmentation 
      join segmentation_feature using (segmentation_id) 
    group by 
      segmentation_id, 
      assignment
  ;
  ~
);

$funcgen_adaptor->dbc->do(
  q~
  insert into segmentation_statistic 
    (
      segmentation_id,
      label,
      statistic, 
      value
    )
    select 
      null,
      assignment,
      'num_states',
      count(distinct state)
    from 
      segmentation 
      join segmentation_feature using (segmentation_id) 
    group by 
      assignment
  ;
  ~
);

$funcgen_adaptor->dbc->do(
  q~
  insert into segmentation_statistic 
    (
      segmentation_id,
      label,
      statistic, 
      value
    )
    select 
      null,
      assignment, 
      'average_length',
      avg(length)
    from 
      segmentation
      join segmentation_feature using (segmentation_id) 
    group by 
      assignment
  ;
  ~
);

# So there are ctcf statistics even for segmentations that don't have ctcf.
$logger->info("Inserting defaults for missing data\n");

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

my @statistics = qw(
  average_length
);

my $segmentation_adaptor = $funcgen_adaptor->get_SegmentationAdaptor;
my $segmentations = $segmentation_adaptor->fetch_all;

my $segmentation_statistic_adaptor 
  = $funcgen_adaptor->get_SegmentationStatisticAdaptor;

foreach my $segmentation (@$segmentations) {
  foreach my $segmentation_assignment (@segmentation_assignments) {
    STATISTIC:
    foreach my $statistic (@statistics) {
    
      $logger->info("Checking $statistic, " . $segmentation->name . ", $segmentation_assignment\n");

      my $segmentation_statistic
        = $segmentation_statistic_adaptor
          ->_statistic_exists(
            $statistic, 
            $segmentation, 
            $segmentation_assignment
          );
      if (@$segmentation_statistic > 0) {
        next STATISTIC;
      }
      my $segmentation_id = $segmentation->dbID;
      my $sql = 
        qq~
        insert into segmentation_statistic 
          (
            segmentation_id,
            label,
            statistic, 
            value
          ) values (
            $segmentation_id,
            '$segmentation_assignment',
            '$statistic',
            0
          )
        ;
        ~;
      $logger->info("$sql\n");
      $funcgen_adaptor->dbc->do($sql);
    }
  }
}

$logger->finish_log;
exit(0);
