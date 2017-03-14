#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Getopt::Long;

=head1 import_parse_exonerate.pl

perl scripts/array_mapping/import_parse_exonerate.pl \
  --exonerate_file /nfs/nobackup/ensembl/mnuhn/array_mapping/temp/homo_sapiens/probe_chunks/probe_chunk_2.fasta_genomic.exonerate.txt \
  --max_allowed_mismatches_per_hit 0 \
  > /nfs/nobackup/ensembl/mnuhn/array_mapping/temp/homo_sapiens/probe_chunks/probe_chunk_2.fasta_genomic.exonerate_parsed.txt

=cut

my $exonerate_file;
my $parsed_output;
my $max_allowed_mismatches_per_hit;

GetOptions (
   'exonerate_file=s'                 => \$exonerate_file,
   'parsed_output=s'                  => \$parsed_output,
   'max_allowed_mismatches_per_hit=s' => \$max_allowed_mismatches_per_hit,
);

if (! -e $exonerate_file) {
  die("Can't find probe file ${exonerate_file}!");
}

open my $exonerate_fh, '<', $exonerate_file;

my $current_probe_seq_id;
my $current_list_of_hits_for_this_probe;

# This makes debugging easier, if someone has to inspect the output of this 
# script.
#
$Data::Dumper::Sortkeys = 1;

HIT: while (my $current_line = <$exonerate_fh>) {

  next unless ($current_line=~/^RESULT:/);
  chomp;

  my @f = split ' ', $current_line;
  
  my %hit;
  
  $hit{probe_seq_id}   = $f[1];
  $hit{q_start}        = $f[2];
  $hit{q_end}          = $f[3];
  $hit{q_strand}       = $f[4];
  $hit{t_id}           = $f[5];
  $hit{t_start}        = $f[6];
  $hit{t_end}          = $f[7];
  $hit{t_strand}       = $f[8];
  $hit{score}          = $f[9];
  $hit{perc_id}        = $f[10];
  $hit{q_length}       = $f[11];
  $hit{t_length}       = $f[12];
  $hit{mismatch_count} = $f[13];
  $hit{scores}         = $f[14];
  
  my $cigar_line  = compute_cigar_line(\%hit);
  my $hit_ensembl = compute_ensembl_coordinates(\%hit);
  
  $hit_ensembl->{cigar_line}       = $cigar_line;
  $hit_ensembl->{total_mismatches} = compute_total_mismatches(\%hit);

  next HIT if ($hit_ensembl->{total_mismatches} > $max_allowed_mismatches_per_hit);
  
  if ($current_probe_seq_id != $hit{probe_seq_id}) {
  
    # This will be empty in the first iteration and any probe that doesn't 
    # make any hits.
    #
    if ($current_list_of_hits_for_this_probe) {
      print Dumper($current_list_of_hits_for_this_probe);
    }
    $current_list_of_hits_for_this_probe = [];
    $current_probe_seq_id = $hit_ensembl->{probe_seq_id};
    
  }
  push @$current_list_of_hits_for_this_probe, $hit_ensembl;
}

sub compute_match_length {
  my $hit = shift;  
  return $hit->{q_end} - $hit->{q_start};
}

sub compute_align_mismatch {
  my $hit = shift;  
  return $hit->{q_length} - compute_match_length($hit);
}

sub compute_total_mismatches {
  my $hit = shift;  
  return $hit->{mismatch_count} + compute_align_mismatch($hit);
}

sub compute_ensembl_coordinates {
  my $original_hit = shift;

  my @soft_cigar_line;
  
  # Make a copy to avoid messing up the original.
  #
  my %tmp = %$original_hit;
  my $hit = \%tmp;

  my $q_length = $hit->{q_length};
  my $q_start  = $hit->{q_start};
  my $q_end    = $hit->{q_end};
  my $t_start  = $hit->{t_start};
  my $t_end    = $hit->{t_end};
  my $t_strand = $hit->{t_strand};

  # Convert to Ensembl coordinates
  if($hit->{t_strand} eq '+') {
    $hit->{t_start} += 1;
  } else {
    $hit->{t_end}   += 1;
  }

  my $align_mismatch = compute_align_mismatch($hit);
  
  my $query_does_not_match_at_beginning = $align_mismatch && ($q_length == $q_end);
  
  if($query_does_not_match_at_beginning) {
    $t_start = ($t_strand eq '+') ? ($t_start - $q_start) : ($t_start + $q_start);
    $hit->{t_start} = $t_start;
  }
  if($align_mismatch != $q_start) {
    #Add end mismatch if
    #not accounted for by 5' mismatch
    my $three_mismatch = ($q_length - $q_end);
    push @soft_cigar_line, $three_mismatch.'X';
    $t_end = ($t_strand eq '+') ? ($t_end + $three_mismatch) : ($t_end - $three_mismatch);
    $hit->{t_end} = $t_end;
  }

  if (
    ($hit->{q_strand} eq '+') && ($hit->{t_strand} eq '-')
  ) {
    ($hit->{t_end}, $hit->{t_start}) = ($hit->{t_start}, $hit->{t_end});
  }
  eval {
    $hit->{t_strand} = exonerate_strand_to_ensembl_strand($hit->{t_strand});
    # Should always be +
    $hit->{q_strand} = exonerate_strand_to_ensembl_strand($hit->{q_strand});
  };
  if ($@) {
    die(Dumper($hit));
  }
  return $hit;
}

sub exonerate_strand_to_ensembl_strand {

  my $exonerate_strand = shift;
  
  if ($exonerate_strand eq '+') {
    return 1
  }
  if ($exonerate_strand eq '-') {
    return -1
  }
  die("Unrecognised target strand symbol: $exonerate_strand");
}

sub compute_cigar_line {

  my $original_hit = shift;

  my @soft_cigar_line;
  
  # Make a copy to avoid messing up the original.
  #
  my %tmp = %$original_hit;
  my $hit = \%tmp;

  my $q_length = $hit->{q_length};
  my $q_start  = $hit->{q_start};
  my $q_end    = $hit->{q_end};
  my $q_strand = $hit->{q_strand};
  my $t_start  = $hit->{t_start};
  my $t_end    = $hit->{t_end};
  my $t_strand = $hit->{t_strand};

  my $align_mismatch = compute_align_mismatch($hit);
  
  my $query_does_not_match_at_beginning = $align_mismatch && ($q_length == $q_end);
  
  if($query_does_not_match_at_beginning) {
  
    push @soft_cigar_line, $q_start.'X' if $q_start;#set this to the value of start if not 0
  }

  my $total_mismatches = compute_total_mismatches($hit);
  
  if($total_mismatches) {

    # Looks like this:
    #
    # 'scores' => 'scores:0:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:0:',
    #
    my @scores = split ':', $hit->{scores};
    
    # Removes "scores"
    shift @scores;
    # Removes "0"
    shift @scores;
    # Removes last 0
    pop @scores;

    my $number_of_identical_adjacent_scores = 0;
    
    # Always starts with a match otherwise it wouldn't have been reported.
    my $prev_score = 5;

    foreach my $tscore (@scores) {
    
      if($tscore == $prev_score) {
        $number_of_identical_adjacent_scores += 1;
      } else {
        my $tmp = ($prev_score == 5) ?  $number_of_identical_adjacent_scores.'=' :  $number_of_identical_adjacent_scores.'X';
        push @soft_cigar_line, $tmp;
        $number_of_identical_adjacent_scores = 1;
        $prev_score       = $tscore;
      }
    }
    #handle last align length
    my $tmp = ($prev_score == 5) ? $number_of_identical_adjacent_scores.'=' : $number_of_identical_adjacent_scores.'X';
    push @soft_cigar_line, $tmp;
  }
  
  if ( ($q_strand eq '+') && ($t_strand eq '-') ) {
    @soft_cigar_line = reverse(@soft_cigar_line);
  }

  my $cigar_string = join(':', @soft_cigar_line) || compute_match_length($hit) . '=';
  return $cigar_string;
}

