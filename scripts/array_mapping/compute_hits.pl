#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Getopt::Long;
use Hash::Util qw( lock_hash lock_keys );

=head1

compute_hits.pl \
  --probeset_sizes_file   /nfs/nobackup/ensembl/mnuhn/probe2transcript/probeset_sizes.pl   \
  --transcript_info_file  /nfs/nobackup/ensembl/mnuhn/probe2transcript/transcript_info.pl \
  --probeset_transcript_assignments /nfs/nobackup/ensembl/mnuhn/probe2transcript/probeset_transcript_assignments.pl \
  --probeset_transcript_rejections  /nfs/nobackup/ensembl/mnuhn/probe2transcript/probeset_transcript_rejections.pl \
  --probe_transcript_assignments    /nfs/nobackup/ensembl/mnuhn/probe2transcript/probe_transcript_assignments.pl \

=cut

# Constants:

my $mapping_threshold = 0.5;

my $transcript_info_file;
my $probeset_sizes_file;
my $probeset_transcript_assignments_file;
my $probeset_transcript_rejections_file;
my $probe_transcript_assignments_file;

GetOptions (
   'transcript_info_file=s'                 => \$transcript_info_file,
   'probeset_sizes_file=s'                  => \$probeset_sizes_file,
   'probeset_transcript_assignments_file=s' => \$probeset_transcript_assignments_file,
   'probeset_transcript_rejections_file=s'  => \$probeset_transcript_rejections_file,
   'probe_transcript_assignments_file=s'    => \$probe_transcript_assignments_file,
);

use Bio::EnsEMBL::Utils::Logger;
my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

use Bio::EnsEMBL::Funcgen::Parsers::DataDumper;
my $parser = Bio::EnsEMBL::Funcgen::Parsers::DataDumper->new;

my $probeset_sizes = [];

$logger->info("Reading probeset_sizes from " . $probeset_sizes_file . ".\n");

$parser->parse({
  data_dumper_file => $probeset_sizes_file,
  call_back        => sub {
    my $probeset_sizes_for_one_array = shift;
    push @$probeset_sizes, $probeset_sizes_for_one_array;
  },
});

$logger->info("Got probeset_sizes from " . @$probeset_sizes . " arrays.\n");

$logger->info("Writing transcript assignments to " . $probeset_transcript_assignments_file . ".\n");
$logger->info("Writing transcript rejections to "  . $probeset_transcript_rejections_file  . ".\n");

open my $probeset_transcript_assignments_fh, '>', $probeset_transcript_assignments_file;
open my $probeset_transcript_rejections_fh,  '>', $probeset_transcript_rejections_file;
open my $probe_transcript_assignments_fh,    '>', $probe_transcript_assignments_file;

my $compute_hits_parameters = {
  probeset_transcript_assignments_fh => $probeset_transcript_assignments_fh,
  probeset_transcript_rejections_fh  => $probeset_transcript_rejections_fh,
  probe_transcript_assignments_fh    => $probe_transcript_assignments_fh,
  probeset_sizes                     => $probeset_sizes,
};

$parser->parse({
  data_dumper_file => $transcript_info_file,
  call_back        => sub {
    my $transcript_info = shift;
    $compute_hits_parameters->{transcript_info} = $transcript_info;
    assign_probes_to_transcripts($compute_hits_parameters);
  },
});

$probeset_transcript_assignments_fh->close;
$probeset_transcript_rejections_fh->close;
$probe_transcript_assignments_fh->close;

$logger->finish_log;

sub assign_probes_to_transcripts {
  
  my $param = shift;
  lock_keys(%$param);

  my $transcript_info                    = $param->{transcript_info};
  my $all_probeset_sizes                 = $param->{probeset_sizes};
  my $probeset_transcript_assignments_fh = $param->{probeset_transcript_assignments_fh};
  my $probeset_transcript_rejections_fh  = $param->{probeset_transcript_rejections_fh};
  my $probe_transcript_assignments_fh    = $param->{probe_transcript_assignments_fh};
  
  lock_hash(%$transcript_info);
  
  my $transcript_stable_id = $transcript_info->{transcript_stable_id};
  my $probe_hits_by_array  = $transcript_info->{probe_hits_by_array};
  
  my @all_array_names = keys %$probe_hits_by_array;
  
  ARRAY:
  foreach my $current_array_name (@all_array_names) {
  
    next ARRAY if (! exists $probe_hits_by_array->{$current_array_name});
  
    my $hit_summary = $probe_hits_by_array->{$current_array_name};
    
#     print Dumper($hit_summary);
#     
#     if (ref $hit_summary ne 'HASH') {
#       print "Not a hash reference for $current_array_name\n";
#     } else {
#       print "Got a " . (ref $hit_summary) . " for $current_array_name\n";
#     }
    
    lock_hash(%$hit_summary);
    
    my $is_probeset_array = $hit_summary->{probeset_array};
    
    my $probeset_sizes = undef;
    
    if ($is_probeset_array) {
      ($probeset_sizes) = grep { $_->{array_name} eq $current_array_name } @$all_probeset_sizes;
      if (! defined $probeset_sizes) {
        die;
      }
      use Hash::Util qw( lock_keys );
      lock_keys(%$probeset_sizes);
      (
        my $assignments, 
        my $rejections
      ) = assign_probsets_with_sufficient_probe_matches({
          transcript_stable_id => $transcript_stable_id,
          hit_summary          => $hit_summary,
          probeset_sizes       => $probeset_sizes->{probeset_sizes},
      });
      write_results_as_tsv($probeset_transcript_assignments_fh, $assignments);
      write_results_as_tsv($probeset_transcript_rejections_fh,  $rejections);
    } else {
        my $assignments = assign_all_probes({
          transcript_stable_id => $transcript_stable_id,
          hit_summary          => $hit_summary,
      });
      write_results_as_tsv($probe_transcript_assignments_fh,$assignments);
    }
  }
  return;
}

sub write_results_as_tsv {
  my $fh   = shift;
  my $hash = shift;
  foreach my $key (keys %$hash) {
    $fh->print($key . "\t" . $hash->{$key} . "\n");
  }
}

sub assign_probsets_with_sufficient_probe_matches {

  my $param = shift;
  lock_hash(%$param);

  my $transcript_stable_id = $param->{transcript_stable_id};
  my $hit_summary          = $param->{hit_summary};
  my $probeset_sizes       = $param->{probeset_sizes};

  # Count hits between object and transcript
  my $hits;
  my $hits_are_sufficient;
  my %assignments;
  my %rejections;
  
  my $probesets_matching_transcript = $hit_summary->{probesets};

  #  probesets_matching_transcript is a hash with each entry looking like this:
  # 
  #   '755018' => {
  #     'num_mismatches' => 0,
  #     'num_probe_features' => 21,
  #     'probe_id' => {
  #                     '6573852' => {
  #                                     'count' => 5
  #                                   },
  #                     '6573850' => {
  #                                     'count' => 5
  #                                   },
  #                     '6573851' => {
  #                                     'count' => 5
  #                                   },
  #                     '6573849' => {
  #                                     'count' => 6
  #                                   }
  #                   }
  #  },
  
  my @probeset_ids = keys %$probesets_matching_transcript;
  
  foreach my $current_probeset_id (@probeset_ids) {
  
    # The number of probes from that probeset that have a probefeature
    # on this transcript.
    #
    my $hits = scalar keys %{$probesets_matching_transcript->{$current_probeset_id}};
    
    my $probeset_size = $probeset_sizes->{$current_probeset_id};
    
    my $hit_is_sufficient = ($hits / $probeset_size) >= $mapping_threshold;
    
    if ($hit_is_sufficient) {
      $assignments{$current_probeset_id} = $transcript_stable_id;
    } else {
      $rejections{$current_probeset_id} = $transcript_stable_id;
    }
  }
  return (\%assignments, \%rejections);
}

sub assign_all_probes {

  my $param = shift;
  lock_hash(%$param);

  my $transcript_stable_id = $param->{transcript_stable_id};
  my $hit_summary          = $param->{hit_summary};

  my %assignments;
  
  # The hash looks like this;
  #
  # 'OneArray' => {
  #   'array_name' => 'OneArray',
  #   'probe' => {
  #       '13064' => [
  #                     {
  #                       'has_mismatches' => 0,
  #                       'annotation' => 'exon'
  #                     },
  #                     {
  #                       'has_mismatches' => 0,
  #                       'annotation' => 'exon'
  #                     },
  #                     {
  #                       'has_mismatches' => 0,
  #                       'annotation' => 'exon'
  #                     }
  #                   ],
  #       '25792' => [
  #                     {
  #                       'has_mismatches' => 0,
  #                       'annotation' => 'exon'
  #                     },
  #                     {
  #                       'has_mismatches' => 0,
  #                       'annotation' => 'exon'
  #                     },
  #                     {
  #                       'has_mismatches' => 0,
  #                       'annotation' => 'exon'
  #                     },
  #                     {
  #                       'has_mismatches' => 0,
  #                       'annotation' => 'exon'
  #                     },
  #                     {
  #                       'has_mismatches' => 0,
  #                       'annotation' => 'exon'
  #                     }
  #                   ]
  #     },
  #   'probeset_array' => undef
  # },

  my @probe_ids_with_valid_hits = keys %{$hit_summary->{probe}};
  
  foreach my $current_probe_id (@probe_ids_with_valid_hits) {
    $assignments{$current_probe_id} = $transcript_stable_id;
  }
  return \%assignments;
}
