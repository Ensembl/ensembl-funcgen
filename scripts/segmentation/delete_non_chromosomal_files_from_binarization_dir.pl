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

  delete_non_chromosomal_files_from_binarization_dir.pl \
    --species  homo_sapiens \
    --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb7.pm \
    --directory /hps/nobackup/production/sds-flicek-bp/blueprint_fastq_files/mnuhn_regbuild_pipeline/rb_human_encode_blueprint_and_tfs/temp_dir/segmentation/homo_sapiens/binarization/encode/no_ctcf/

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
    "directory|d=s",
 );

use Hash::Util qw( lock_keys );
lock_keys( %options );

my $species   = $options{'species'};
my $registry  = $options{'registry'};
my $directory = $options{'directory'};

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

$logger->info("registry  = " . $registry  . "\n");
$logger->info("species   = " . $species   . "\n");
$logger->info("directory = " . $directory . "\n");

Bio::EnsEMBL::Registry->load_all($registry);
my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor( $species, 'core', 'slice' );
if (! defined $slice_adaptor) {
  die("Can't get coord system adaptor! Please configure your registry accordingly.")
}

opendir(my $dh, $directory) || die "Can't opendir $directory: $!";
my @files = grep { -f $_ } map { "$directory/$_" } readdir($dh);
closedir $dh;

foreach my $current_file (@files) {
  #$logger->info("Checking $current_file\n");
  
  my $is_on_chromosome_level = check_file($current_file);
  
  if ($is_on_chromosome_level) {
    $logger->info("Keep $current_file\n");
  } else {
    $logger->info("Deleting $current_file\n");
    unlink($current_file)
  }

}

$logger->finish_log;
exit(0);

sub check_file {

  my $file = shift;

  open my $fh, '<', $file || die("Can't open $file!");

  my $first_line = <$fh>;
  chomp $first_line;
  
  (
    my $epigenome_name,
    my $seq_region_name,
  )
    = split "\t", $first_line;

  my $slice = $slice_adaptor->fetch_by_region(undef , $seq_region_name);
  if (! defined $slice) {
    die("Can't find sequence region $seq_region_name in database!");
  }
  my $is_on_chromosome_level = $slice->coord_system->name eq 'chromosome';
  $fh->close;
  return $is_on_chromosome_level;
}

