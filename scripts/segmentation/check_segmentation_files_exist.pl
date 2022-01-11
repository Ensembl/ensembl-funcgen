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

  check_segmentation_files_exist.pl \
    --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb7.pm \
    --species homo_sapiens \
    --db_file_path /hps/nobackup/production/sds-flicek-bp/blueprint_fastq_files/mnuhn_regbuild_pipeline/rb_human_merged_old_and_new/dbfiles/

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
    "db_file_path|d=s",
 );

use Hash::Util qw( lock_keys );
lock_keys( %options );

my $species      = $options{'species'};
my $registry     = $options{'registry'};
my $db_file_path = $options{'db_file_path'};

Bio::EnsEMBL::Registry->load_all($registry);
my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');

my $coord_system_adaptor = Bio::EnsEMBL::Registry->get_adaptor( $species, 'Core', 'CoordSystem' );
if (!$coord_system_adaptor) {
  die("Can't get coord system adaptor! Please configure your registry accordingly.")
}
my ($cs) = @{$coord_system_adaptor->fetch_all()};
my $assembly = $cs->version();
if (!$assembly) {
  die("Can't work out assembly for $species!")
}

my $segmentation_file_adaptor = $funcgen_adaptor->get_SegmentationFileAdaptor;
my $segmentation_files = $segmentation_file_adaptor->fetch_all;

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

SEGMENTATION_FILE:
foreach my $segmentation_file (@$segmentation_files) {

  my $full_path = join '/', $db_file_path, $species, $assembly, $segmentation_file->file;
  
  if (-e $full_path) {
    $logger->info("ok - file exists: $full_path\n" );
    next SEGMENTATION_FILE;
  }
  $logger->error("not ok - file does not exist: $full_path\n" );
  exit(1);
}

$logger->finish_log;
exit(0);


