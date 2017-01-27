#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Getopt::Long;
use Hash::Util qw( lock_hash lock_keys );

=head1

time create_probe_to_transcript_descriptions.pl \
  --array_name SurePrint_GPL16709_4x44k \
  --probe_transcript_hits_by_array_file /nfs/nobackup/ensembl/mnuhn/array_mapping/temp/oryctolagus_cuniculus/probe2transcript/probe_transcript_hits_file.pl \
  --probe_to_transcript_file         /nfs/nobackup/ensembl/mnuhn/array_mapping/temp/oryctolagus_cuniculus/probe2transcript/probe_to_transcript_file

=cut

my $array_name;
my $probe_transcript_hits_by_array_file;
my $probe_to_transcript_file;

GetOptions (
   'array_name=s'                          => \$array_name,
   'probe_transcript_hits_by_array_file=s' => \$probe_transcript_hits_by_array_file,
   'probe_to_transcript_file=s'            => \$probe_to_transcript_file,
);

# my $w = ['one', 'two', 'three', 'four'];
# my $x = ['one', 'two', 'three'];
# my $y = ['one', 'two'];
# my $z = ['one'];
# 
# print hit_site_to_english($w);
# print "\n";
# print hit_site_to_english($x);
# print "\n";
# print hit_site_to_english($y);
# print "\n";
# print hit_site_to_english($z);
# print "\n";
# 
# exit;

open my $probe_to_transcript_fh, '>', $probe_to_transcript_file;

use Bio::EnsEMBL::Utils::Logger;
my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

use Bio::EnsEMBL::Funcgen::Parsers::DataDumper;
my $parser = Bio::EnsEMBL::Funcgen::Parsers::DataDumper->new;

$logger->info("Reading probe transcript hits file " . $probe_transcript_hits_by_array_file . ".\n");

# my $probe_transcript_hits_by_array;

$parser->parse({
  data_dumper_file => $probe_transcript_hits_by_array_file,
  call_back        => sub {
#     my $x = shift;
#     if ($probe_transcript_hits_by_array) {
#       # There should only be one hash in this.
#       die;
#     }
#     $probe_transcript_hits_by_array = $x;
    my $probe_transcript_hits_by_array = shift;
    create_probe_transcript_description(
      $probe_transcript_hits_by_array, 
      $probe_to_transcript_fh
    );
  },
});

$logger->info("Done reading probe transcript hits file\n");

$probe_to_transcript_fh->close;

$logger->finish_log;
$logger->info("All done.\n");


# Iterate over something that has this structure:
# 
# $VAR1 = {
#           'SurePrint_GPL16709_4x44k' => {
#                                           '8563' => {
#                                                       'ENSOCUT00000013244' => [
#                                                                                 {
#                                                                                   'has_mismatches' => 0,
#                                                                                   'annotation' => 'exon'
#                                                                                 },
#                                                                                 {
#                                                                                   'has_mismatches' => 0,
#                                                                                   'annotation' => 'exon'
#                                                                                 }
#                                                                               ]
#                                                     },
#                                           '28336' => {
#                                                        'ENSOCUT00000012968' => [
#                                                                                  {
#                                                                                    'has_mismatches' => 0,
#                                                                                    'annotation' => '3\' flank'
#                                                                                  }
#                                                                                ]
#                                                      },
#                                           '8434' => {
#                                                       'ENSOCUT00000015991' => [
#                                                                                 {
#                                                                                   'has_mismatches' => 0,
#                                                                                   'annotation' => 'exon/5\' flank boundary'
#                                                                                 }
#                                                                               ]
#                                                     },
sub create_probe_transcript_description {
  
  my $probe_transcript_hits_by_array = shift;
  my $probe_to_transcript_fh         = shift;
  
  my @array_names = keys %$probe_transcript_hits_by_array;

  foreach my $current_array (@array_names) {

    my $probe_transcript_hits = $probe_transcript_hits_by_array->{$current_array};
    my @probe_names = keys %$probe_transcript_hits;
    
    foreach my $current_probe_name (@probe_names) {
    
      my @final_probe_assignments;
    
      my $current_probe_hits = $probe_transcript_hits->{$current_probe_name};
      my @stable_ids_mapped_to = keys %$current_probe_hits;
      
      foreach my $current_stable_id (@stable_ids_mapped_to) {
        
        my $match_description = $current_probe_hits->{$current_stable_id};
        my $num_probe_features_mapped = scalar @{$current_probe_hits->{$current_stable_id}};
        
        push @final_probe_assignments, {
          num_probe_features_mapped => $num_probe_features_mapped,
          stable_id                 => $current_stable_id,
          match_description         => $match_description,
        };
      }

      my $num_stable_ids_assigned_minus_one = -1 + scalar @final_probe_assignments;
      my $other_transcripts_text;
      if ($num_stable_ids_assigned_minus_one==0) {
        $other_transcripts_text = "Matches uniquely to this transcript.";
      }
      if ($num_stable_ids_assigned_minus_one==1) {
        $other_transcripts_text = "Matches one other transcript.";
      }
      if ($num_stable_ids_assigned_minus_one>1) {
        $other_transcripts_text = "Matches $num_stable_ids_assigned_minus_one other transcripts.";
      }

      foreach my $current_final_probe_assignment (@final_probe_assignments) {
      
        my $num_probe_features_mapped = $current_final_probe_assignment->{num_probe_features_mapped};
        my $stable_id                 = $current_final_probe_assignment->{stable_id};
        my $match_description         = $current_final_probe_assignment->{match_description};
        
        use List::Util qw( uniqstr );
        my @hit_site = uniqstr sort map { $_->{annotation} } @$match_description;
        
        my $hit_site_enumeration = hit_site_to_english(\@hit_site);
      
        my $probe_transcript_description = "Matches ${hit_site_enumeration}. $other_transcripts_text";
        
        $probe_to_transcript_fh->print(
          join "\t", $current_probe_name, $stable_id, $probe_transcript_description
        );
        $probe_to_transcript_fh->print("\n");
      }
    }
  }
}

sub hit_site_to_english {

  my $hit_site = shift;
  
  if (@$hit_site == 1) {
    return $hit_site->[0];
  }
  if (@$hit_site == 2) {
    my $first  = shift @$hit_site;
    my $second = shift @$hit_site;
    return "$first and $second";
  }
  if (@$hit_site > 2) {
    my $first = shift @$hit_site;
    my @rest  = @$hit_site;
    return $first . ', ' . hit_site_to_english(\@rest);
  }
  die "Got empty list!";
}


