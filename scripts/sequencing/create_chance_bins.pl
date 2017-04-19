#!/usr/bin/env perl
=head1

=head2 Description

  Script for generating a bin file that is required to run the quality
  check CHANCE for ChIP-seq data.
  
  This script it required by the ERSA pipeline.

=head2 Example

perl test_create_interval_list.pl \
  --species mus_musculus \
  --epigenome_gender hermaphrodite \
  --assembly GRCm38 \
  --outputfile deleteme.bed \
  --reference_data_root_dir /lustre/scratch109/ensembl/funcgen/refbuilder

=cut
use strict;
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::Hive::RefBuildFileLocator;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Logger;

my %options;
GetOptions (
    \%options,
    "species|s=s",
    "epigenome_gender|e=s",
    "assembly|a=s",
    "outputfile|o=s",
    "reference_data_root_dir|r=s",
 );

use Hash::Util qw( lock_keys );
lock_keys( %options );

my $outputfile              = $options{'outputfile'};
my $species                 = $options{'species'};
my $epigenome_gender        = $options{'epigenome_gender'};
my $assembly                = $options{'assembly'};
my $reference_data_root_dir = $options{'reference_data_root_dir'};

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

my $bwa_index_locator = Bio::EnsEMBL::Funcgen::Hive::RefBuildFileLocator->new;

my $chromosome_lengths_relative = $bwa_index_locator->locate({
  species          => $species,
  epigenome_gender => $epigenome_gender,
  assembly         => $assembly,
  file_type        => 'chromosome_lengths_by_species_assembly',
});

my $chromosome_lengths = $reference_data_root_dir . '/' . $chromosome_lengths_relative;

open(my $IN,  '<', $chromosome_lengths) 
  or $logger->error("Couldn't open file $chromosome_lengths:\n$!");
$logger->info("Reading chromosome lengths from $chromosome_lengths\n");

open(my $OUT, '>', $outputfile) or die "$!";
$logger->info("Writing bins to $outputfile\n");

while (my $current_line = <$IN>) {

  chomp $current_line;
  my @s = split "\t", $current_line;
  
  my $seq_region_name   = $s[0];
  my $seq_region_length = $s[1];
  
  my $bin_start = 0;

  for($bin_start =0; $seq_region_length - $bin_start > 1000; $bin_start += 1000) {
    print $OUT "$seq_region_name\t$bin_start\t".($bin_start+1000)."\n";
  }
  # Print last line
  print $OUT "$seq_region_name\t$bin_start\t$seq_region_length\n";
}

$OUT->close;

use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd );
run_system_cmd("bedSort $outputfile $outputfile", undef, 1);

$logger->finish_log;
