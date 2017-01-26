#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Getopt::Long;

=head1
=head2 Description

  This script takes the probe_transcript_hits_file file from 
  compute_probeset_transcript_hits.pl and the probeset_sizes file.
  
=cut

=head1

bsub -q production-rh7 -M64000 -R"select[mem>64000] rusage[mem=64000]" -Is $SHELL

time create_probeset_to_transcript_descriptions.pl \
  --array_name HuGene-2_0-st-v1 \
  --probeset_sizes                      /nfs/nobackup/ensembl/mnuhn/array_mapping/temp/homo_sapiens/probeset_sizes.pl \
  --probeset_transcript_hits_by_array_file /nfs/nobackup/ensembl/mnuhn/array_mapping/temp/homo_sapiens/probe_transcript_hits_file.HuGene-2_0-st-v1.pl \
  --probeset_to_transcript_file         /nfs/nobackup/ensembl/mnuhn/array_mapping/temp/homo_sapiens/probeset_to_transcript_file.pl

0 warnings. Runtime: 0h 2min 52sec [2017-01-25 11:42:27, mem 6.2 Gb]

=cut

# Constants:
my $mapping_threshold = 0.5;

my $array_name;
my $probeset_sizes_file;
my $probeset_transcript_hits_by_array_file;
my $probeset_to_transcript_file;

GetOptions (
   'array_name=s'                 => \$array_name,
   'probeset_sizes_file=s'        => \$probeset_sizes_file,
   'probeset_transcript_hits_by_array_file=s' => \$probeset_transcript_hits_by_array_file,
   'probeset_to_transcript_file=s'      => \$probeset_to_transcript_file,
);

open my $probeset_to_transcript_fh, '>', $probeset_to_transcript_file;

use Bio::EnsEMBL::Utils::Logger;
my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

use Bio::EnsEMBL::Funcgen::Parsers::DataDumper;
my $parser = Bio::EnsEMBL::Funcgen::Parsers::DataDumper->new;

my %probeset_sizes;
$logger->info("Reading probeset_sizes from " . $probeset_sizes_file . ".\n");

$parser->parse({
  data_dumper_file => $probeset_sizes_file,
  call_back        => sub {
    my $probeset_sizes_for_one_array = shift;
    $probeset_sizes{$probeset_sizes_for_one_array->{array_name}} = $probeset_sizes_for_one_array;
  },
});

$logger->info("Got probeset_sizes from " . (scalar (keys %probeset_sizes)) . " arrays.\n");
$logger->info("Reading probe transcript hits file " . $probeset_transcript_hits_by_array_file . ".\n");

my $probe_transcript_hits_by_array;

$parser->parse({
  data_dumper_file => $probeset_transcript_hits_by_array_file,
  call_back        => sub {
    my $x = shift;
    if ($probe_transcript_hits_by_array) {
      # There should only be one hash in this.
      die;
    }
    $probe_transcript_hits_by_array = $x;
  },
});

$logger->info("Done reading probe transcript hits file\n");

my @array_names = keys %$probe_transcript_hits_by_array;

$logger->info("Found " . scalar @array_names . " arrays in the probe transcript hits file.\n");

foreach my $current_array (@array_names) {

  $logger->info("The current array is " . $current_array . " .\n");

  my $probe_transcript_hits = $probe_transcript_hits_by_array->{$current_array};
  my @probeset_names = keys %$probe_transcript_hits;

  $logger->info("It has " . scalar @probeset_names . " probesets.\n");
  
  my $this_arrays_probeset_sizes = $probeset_sizes{$current_array}->{probeset_sizes};
  
  
  foreach my $current_probeset_name (@probeset_names) {
  
    my @final_probeset_assignments;
  
    my $current_probeset_hits = $probe_transcript_hits->{$current_probeset_name};
    my @stable_ids_mapped_to = keys %$current_probeset_hits;
    
    my $current_probeset_size = $this_arrays_probeset_sizes->{$current_probeset_name};
    
    foreach my $current_stable_id (@stable_ids_mapped_to) {
      
      my $num_probes_mapped = keys %{$current_probeset_hits->{$current_stable_id}->{probe_id}};
      
      my $probeset_match_accepted = ($num_probes_mapped / $current_probeset_size) >= $mapping_threshold;
      
      if ($probeset_match_accepted) {
        push @final_probeset_assignments, {
          num_probes_mapped     => $num_probes_mapped,
          current_probeset_size => $current_probeset_size,
          stable_id             => $current_stable_id,
        };
      }
    }

    my $num_stable_ids_assigned_minus_one = -1 + scalar @final_probeset_assignments;
    
    foreach my $current_final_probeset_assignment (@final_probeset_assignments) {
    
      my $num_probes_mapped     = $current_final_probeset_assignment->{num_probes_mapped};
      my $current_probeset_size = $current_final_probeset_assignment->{current_probeset_size};
      my $stable_id             = $current_final_probeset_assignment->{stable_id};
    
      my $probeset_transcript_description = "$num_probes_mapped out of $current_probeset_size probes from this probeset have been mapped to this transcript. The probeset matches $num_stable_ids_assigned_minus_one other transcripts.";
      
      $probeset_to_transcript_fh->print(
        join "\t", $current_probeset_name, $stable_id, $probeset_transcript_description
      );
      $probeset_to_transcript_fh->print("\n");
    }
  }
}

$probeset_to_transcript_fh->close;

$logger->finish_log;
$logger->info("All done.\n");
