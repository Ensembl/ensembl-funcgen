#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2021] EMBL-European Bioinformatics Institute

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

  load_segmentation_files_to_db.pl \
    --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb7.pm \
    --species homo_sapiens \
    --db_file_path /hps/nobackup/production/sds-flicek-bp/blueprint_fastq_files/mnuhn_regbuild_pipeline/rb_human_merged_old_and_new/dbfiles/ \
    --tempdir foobar

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
    "tempdir|t=s",
    "db_file_path|d=s",
 );

use Hash::Util qw( lock_keys );
lock_keys( %options );

my $species      = $options{'species'};
my $registry     = $options{'registry'};
my $tempdir      = $options{'tempdir'};
my $db_file_path = $options{'db_file_path'};

Bio::EnsEMBL::Registry->load_all($registry);
my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');

$funcgen_adaptor->dbc->do(q( drop table if exists segmentation_feature; ));

$funcgen_adaptor->dbc->do(
q~
  create table segmentation_feature (
    segmentation_feature_id   int(28)      unsigned NOT NULL AUTO_INCREMENT,
    segmentation_file_id      int(23)      unsigned not NULL,
    segmentation_id           int(18)      unsigned not NULL,
    epigenome_id              int(22)      unsigned default NULL,
    seq_region_name           varchar(255) default  NULL,
    start                     int(22)      unsigned default NULL,
    end                       int(22)      unsigned default NULL,
    state                     varchar(255) default  NULL,
    assignment                varchar(255) default  NULL,
    length                    int(22)      unsigned default NULL,
    PRIMARY KEY (segmentation_feature_id)
  );
~
);

my $sth_insert_segmentation_feature = $funcgen_adaptor->dbc->prepare(
  q~
    insert into segmentation_feature (
      segmentation_file_id,
      segmentation_id,
      epigenome_id,
      seq_region_name,
      start,
      end,
      state,
      assignment,
      length
    ) values (
      ?, ?, ?, ?, ?, ?, ?, ?, ?
    );
  ~
);

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

my $file_counter = 0;

SEGMENTATION_FILE:
foreach my $segmentation_file (@$segmentation_files) {

  $file_counter++;

  my $full_path = join '/', $db_file_path, $species, $assembly, $segmentation_file->file;
  
  my $epigenome = $segmentation_file->get_Epigenome;
  
  my $segmentation = $segmentation_file->get_Segmentation;
  
  my $current_tempdir = join '/', 
    $tempdir, 
    $epigenome->production_name,
    $segmentation_file->dbID
  ;
  
  use File::Path qw( make_path );
  make_path($current_tempdir);

  $logger->info("\n");
  $logger->info("This is file number $file_counter\n");
  $logger->info("  Copying $full_path to $current_tempdir\n");
  
  use File::Copy;
  copy($full_path, $current_tempdir) or die("Couldn't copy $full_path to $current_tempdir!");
  
  use File::Basename;
  my $segmentation_file_bed_basename = basename($segmentation_file->file, '.bb') . '.bed';
  
  my $segmentation_file_bed_file_name = "$current_tempdir/$segmentation_file_bed_basename";
  
  my $cmd = "bigBedToBed $full_path $segmentation_file_bed_file_name";
  
  $logger->info("Converting to bed\n");
  $logger->info("    $cmd\n");
  
  my $return_value;

  use Capture::Tiny;
  # Capture:Tiny has weird behavior if 'require'd instead of 'use'd
  # see, for example,http://www.perlmonks.org/?node_id=870439 
  my $stderr = Capture::Tiny::tee_stderr(sub {
      $return_value = system($cmd);
  });

  if ($return_value) {
    use Carp;
    confess(
      "Error running command \n\n" 
      . $cmd . "\n\n"
      . " the error was:\n\n" 
      . $stderr . "\n\n"
    );
  }
  if (! -e $segmentation_file_bed_file_name) {
    confess("$segmentation_file_bed_file_name couldn't be created!");
  }
  
  open my $fh, '<', $segmentation_file_bed_file_name || confess("Couldn't open $segmentation_file_bed_file_name!");
  
  $logger->info("  Loading file\n");
  while (my $bed_line = <$fh>) {
  
    (
      my $seq_region_name,
      my $start,
      my $end,
      my $feature_name,
    ) = split "\t", $bed_line;
    
    (
      my $segmentation_state,
      my $segmentation_label,
      my $no_idea_what_this_number_is
    ) = split '_', $feature_name;
    
    # UCSC coordinates, +1 is not needed.
    my $length = $end - $start;
    
    $segmentation_state =~ s/^E//;
    
    $sth_insert_segmentation_feature->execute(
      $segmentation_file->dbID,
      $segmentation->dbID,
      $epigenome->dbID,
      $seq_region_name,
      $start,
      $end,
      $segmentation_state,
      $segmentation_label,
      $length
    );
  }
  $logger->info("  Done loading file\n");
  close($fh);
}

$logger->info("Done loading all files\n");
$logger->finish_log;
exit(0);
