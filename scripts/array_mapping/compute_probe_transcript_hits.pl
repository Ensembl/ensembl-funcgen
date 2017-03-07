#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Getopt::Long;
use Hash::Util qw( lock_hash lock_keys );

=head1

compute_probe_transcript_hits.pl \
  --array_name                 SurePrint_GPL16709_4x44k \
  --transcript_info_file       /nfs/nobackup/ensembl/mnuhn/array_mapping/temp/oryctolagus_cuniculus/probe2transcript/transcript_info.pl \
  --probe_transcript_hits_file /nfs/nobackup/ensembl/mnuhn/array_mapping/temp/oryctolagus_cuniculus/probe2transcript/probe_transcript_hits_file.pl

=cut

my $array_name;
my $transcript_info_file;
my $probe_transcript_hits_file;

GetOptions (
   'array_name=s'                 => \$array_name,
   'transcript_info_file=s'       => \$transcript_info_file,
   'probe_transcript_hits_file=s' => \$probe_transcript_hits_file,
);

use Bio::EnsEMBL::Utils::Logger;
my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

use Bio::EnsEMBL::Funcgen::Parsers::DataDumper;
my $parser = Bio::EnsEMBL::Funcgen::Parsers::DataDumper->new;

$logger->info("Reading transcript info " . $transcript_info_file . ".\n");

my $num_transcripts_seen        = 0;
my $print_life_sign_after_every = 1000;

open my $probe_transcript_hits_file_fh,    '>', $probe_transcript_hits_file;

my %probe_transcript_hits;

$parser->parse({
  data_dumper_file => $transcript_info_file,
  call_back        => sub {
    my $transcript_info = shift;
    process_transcript_info(
      $transcript_info
    );
  },
});
$logger->info("Done reading transcript info\n");

$logger->info("Writing probe transcript hits from array $array_name\n");

my $probe_transcript_hits_for_current_array = $probe_transcript_hits{$array_name};

my @probe_ids = keys %$probe_transcript_hits_for_current_array;
foreach my $current_probe_id (@probe_ids) {

  my $probe_transcript_hits_for_current_probe = $probe_transcript_hits_for_current_array->{$current_probe_id};

  $probe_transcript_hits_file_fh->print(
    Dumper(
      {
        $array_name => 
          {
            $current_probe_id => $probe_transcript_hits_for_current_probe
          }
      }
    )
  );
}

$logger->info("Done.\n");

$probe_transcript_hits_file_fh->close;

$logger->finish_log;
$logger->info("All done.\n");

sub process_transcript_info {

  my $transcript_info = shift;
  
  $num_transcripts_seen++;
  my $transcript_stable_id = $transcript_info->{transcript_stable_id};
  
  if ($num_transcripts_seen % $print_life_sign_after_every == 0) {
    $logger->info("Seen $num_transcripts_seen transcripts.\n");
  }
  
  my $probe_hits_by_array = $transcript_info->{probe_hits_by_array};
  my @array_names = keys %$probe_hits_by_array;
  
  ARRAY:
  foreach my $current_array_name (@array_names) {
  
    next ARRAY unless($current_array_name eq $array_name);
  
    my $current_probe_hits = $probe_hits_by_array->{$current_array_name};
    my $is_probeset_array = $current_probe_hits->{probeset_array};
    
    if ($is_probeset_array) {
    
      my @probe_set_ids = keys %{$current_probe_hits->{probesets}};
      foreach my $current_probe_set_id (@probe_set_ids) {
      
        my @probe_ids = keys %{$current_probe_hits->{probesets}->{$current_probe_set_id}->{probe_id}};
        
        foreach my $current_probe_id (@probe_ids) {
          my $current_hits = $current_probe_hits->{probesets}->{$current_probe_set_id}->{probe_id}->{$current_probe_id};
          $probe_transcript_hits{$current_array_name}{$current_probe_id}{$transcript_stable_id} = $current_hits;
        }
      }
    } else {
      add_probe_hits_to_probe_transcript_hits(
        $current_array_name,
        $current_probe_hits, 
        $transcript_stable_id, 
        \%probe_transcript_hits
      );
    }
  }
}

sub add_probe_hits_to_probe_transcript_hits {

  my $array_name            = shift;
  my $hits                  = shift;
  my $transcript_stable_id  = shift;
  my $probe_transcript_hits = shift;

  my $probe_hits = $hits->{probe};
  my @probe_ids  = keys %$probe_hits;
  
  foreach my $current_probe_id (@probe_ids) {
    my $current_hits = $probe_hits->{$current_probe_id};
    $probe_transcript_hits->{$array_name}{$current_probe_id}{$transcript_stable_id} = $current_hits;
  }
}

