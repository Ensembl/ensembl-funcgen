#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Getopt::Long;

=head1 project_transcript_hits_to_genome.pl

perl scripts/array_mapping/import_parse_exonerate.pl \
  --exonerate_file /nfs/nobackup/ensembl/mnuhn/array_mapping/temp/homo_sapiens/probe_chunks/probe_chunk_2.fasta_transcript.exonerate.txt \
  --max_allowed_mismatches_per_hit 0 \
  > /nfs/nobackup/ensembl/mnuhn/array_mapping/temp/homo_sapiens/probe_chunks/probe_chunk_2.fasta_transcript.exonerate_parsed.txt

perl scripts/array_mapping/import_create_probe_feature_objects.pl \
  --parsed_data /nfs/nobackup/ensembl/mnuhn/array_mapping/temp/homo_sapiens/probe_chunks/probe_chunk_2.fasta_transcript.exonerate_parsed.txt \
  --analysis_logic_name ProbeAlign_transcript \
  --max_allowed_hits_per_probe 100 \
  > /nfs/nobackup/ensembl/mnuhn/array_mapping/temp/homo_sapiens/probe_chunks/probe_chunk_2.fasta_transcript.probe_features.txt

perl scripts/array_mapping/project_transcript_hits_to_genome.pl \
  --probe_features /nfs/nobackup/ensembl/mnuhn/array_mapping/temp/homo_sapiens/probe_chunks/probe_chunk_2.fasta_transcript.probe_features.txt \
  --registry /nfs/gns/homes/mnuhn/work_dir_probemapping/lib/ensembl-funcgen/registry.pm \
  --species homo_sapiens

perl scripts/array_mapping/project_transcript_hits_to_genome.pl \
  --probe_features /nfs/nobackup/ensembl/mnuhn/array_mapping/temp/homo_sapiens/probe_chunks/5/1/probe_chunk_615.fasta_transcript.probe_features.txt \
  --registry /nfs/gns/homes/mnuhn/work_dir_probemapping/lib/ensembl-funcgen/registry.pm \
  --species homo_sapiens > probe_chunk_615.fasta_transcript.probe_features_projected_from_transcript_coords_nogenedups_.txt

=cut

my $probe_features;
my $registry;
my $species;
my $output_file;

GetOptions (
   'probe_features=s' => \$probe_features,
   'registry=s'       => \$registry,
   'species=s'        => \$species,
   'output_file=s'    => \$output_file,
);

if (! -e $probe_features) {
  die("Can't find probe file ${probe_features}!");
}
if (! $output_file) {
  die("output_file is a mandatory parameter!");
}

open my $output_fh, '>', $output_file || die("Couldn't open file ${output_file}!");

use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->load_all($registry);


my $probe_feature_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'probefeature');
my $analysis_adaptor      = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'analysis');
my $gene_adaptor          = Bio::EnsEMBL::Registry->get_adaptor($species, 'core',    'gene');
my $transcript_adaptor    = Bio::EnsEMBL::Registry->get_adaptor($species, 'core',    'transcript');
my $slice_adaptor         = Bio::EnsEMBL::Registry->get_adaptor($species, 'core',    'slice');

$Data::Dumper::Sortkeys = 1;
# $Data::Dumper::Maxdepth = 3;

my %gene_hits;

my $map_transcript_to_genome = sub {

  my $probe_feature_hash = shift;
  
#   if ($probe_feature_hash->{'q_strand'} != $probe_feature_hash->{'t_strand'}) {  
#     return;
#   }
  
  my $transcript_stable_id = $probe_feature_hash->{t_id};
  
  my $transcript = $transcript_adaptor->fetch_by_stable_id($transcript_stable_id);
  
  if ($probe_feature_hash->{q_start} == 0) {
    $probe_feature_hash->{q_start} = 1;
  }
  
  my $length_before_projecting = $probe_feature_hash->{t_end} - $probe_feature_hash->{t_start};
  
  my $transcript_mapper = Bio::EnsEMBL::TranscriptMapper->new($transcript);
  my @genomic_blocks = $transcript_mapper->cdna2genomic(
    $probe_feature_hash->{t_start},
    $probe_feature_hash->{t_end},
  );
  my $projected_hit = project_hit_to_genomic_coordinates({
    hit            => $probe_feature_hash,
    genomic_blocks => \@genomic_blocks,
    transcript     => $transcript,
  });
  if (! defined $projected_hit) {
    warn("Skipping hit that couldn't be projected:\n" . Dumper($probe_feature_hash));
    return;
  }
  
  my $length_after_projecting = $projected_hit->{t_end} - $projected_hit->{t_start};
  
  my $length_changed_during_projection = $length_before_projecting != $length_after_projecting;
  my $cigar_line_shows_insertion = $projected_hit->{cigar_line} =~ /D/;
  
  if ($length_changed_during_projection && !$cigar_line_shows_insertion) {
    die(
      "Cigar line was not adjusted correctly!"
      . Dumper($projected_hit)
    );
  }

  if ($probe_feature_hash->{'q_strand'} != $probe_feature_hash->{'t_strand'}) {
    $projected_hit->{t_strand} = -1 * $projected_hit->{t_strand}
  }

  # Test, if we have already seen this alignment
  my $gene = $gene_adaptor->fetch_by_transcript_stable_id($transcript_stable_id);
  #The only way of doing this is to test the genomic_start/end and the genomic cigarline with the gene_stable_id and the probe_seq_id
  my $gene_stable_id = $gene->stable_id;
  
  my $genomic_start = $projected_hit->{t_start};
  my $genomic_end   = $projected_hit->{t_end};
  my $cigar_line    = $projected_hit->{cigar_line};
  my $probe_seq_id  = $projected_hit->{probe_seq_id};
  
  if (!$cigar_line_shows_insertion) {
  
    my $transcript = $transcript_adaptor->fetch_by_stable_id($transcript_stable_id);
    
    my $transcript_matched_seq = $transcript->seq;
    
    my $matched_sequence = $transcript_matched_seq->subseq(
      $probe_feature_hash->{t_start},
      $probe_feature_hash->{t_end}
    );

    if ($probe_feature_hash->{t_strand} == -1) {
      use Bio::PrimarySeq;
      $matched_sequence = Bio::PrimarySeq->new( -seq => $matched_sequence )->revcom->seq;
    }
    
    my $probe_sequence_adaptor = Bio::EnsEMBL::Registry->get_adaptor( $species, 'Funcgen', 'ProbeSequence' );
    my $probe_sequence_obj = $probe_sequence_adaptor->fetch_by_dbID($probe_seq_id);
    my $probe_sequence = uc($probe_sequence_obj->sequence);

    my $match_ok = $matched_sequence eq $probe_sequence;

    if (!$match_ok) {
    
      $Data::Dumper::Sortkeys = 1;
      $Data::Dumper::Maxdepth = 3;
      
      die(
        "The probe sequence is not identical to the sequence on the transcript:\n\n"
        . ">" . $transcript->stable_id  . "\n"
        . $transcript->seq->seq . "\n"
        . "\n"
        . "Transcript sequence matched: $matched_sequence\n"
        . "Probe sequence:              $probe_sequence\n\n"
        . Dumper($probe_feature_hash)
        . Dumper($probe_sequence)

      );
    }
#   Commented the check out, because
# 
#   - this is already a bottleneck for the pipeline
#   - It is slowing the analysis down even further.
# 
#     my $probe_feature_passes = check_projection_for_perfect_matches(
#       $probe_seq_id,
#       $genomic_start,
#       $genomic_end,
#       $projected_hit->{t_strand},
#       $projected_hit->{t_id},
#     );
#     
#     if (!$probe_feature_passes) {
#       die(
#         "Probe feature didn't pass the projection test!\n"
#         . "Before:\n"
#         . Dumper($probe_feature_hash)
#         . "After:\n"
#         . Dumper($projected_hit)
#       );
#     }
  }
  my $gene_hit_key = "${gene_stable_id}:${probe_seq_id}:${genomic_start}:${genomic_end}:${cigar_line}";

  if (! exists $gene_hits{$gene_hit_key}) {
    $gene_hits{$gene_hit_key} = undef;
    $output_fh->print(Dumper($projected_hit));
  }
};

sub check_projection_for_perfect_matches {

  my $probe_seq_id             = shift;
  my $probe_feature_start  = shift;
  my $probe_feature_end    = shift;
  my $probe_feature_strand = shift;
  my $transcript_stable_id = shift;

  my $transcript = $transcript_adaptor->fetch_by_stable_id($transcript_stable_id);
  
  if (! defined $transcript) {
    die(
      "Can't find transcript for $transcript_stable_id"
      . Dumper($transcript_adaptor)
    );
  }
  
  my $slice = $transcript->slice;

  my $matched_sequence = $slice->subseq(
    $probe_feature_start,
    $probe_feature_end,
    $probe_feature_strand,
  );
  
  my $probe_sequence_adaptor = Bio::EnsEMBL::Registry->get_adaptor( $species, 'Funcgen', 'ProbeSequence' );
  my $probe_sequence_obj = $probe_sequence_adaptor->fetch_by_dbID($probe_seq_id);
  my $probe_sequence = uc($probe_sequence_obj->sequence);
  
  my $match_ok = $matched_sequence eq $probe_sequence;
  
  my $planned_probe_feature_passed;
  
  if ($match_ok) {
    $planned_probe_feature_passed = 1;
  }
  if (! $match_ok) {
    $planned_probe_feature_passed = undef;
    warn "Not ok: " . $matched_sequence . " != " . $probe_sequence . " " . $probe_feature_strand . "\n";
  }
  return $planned_probe_feature_passed;
}

use Bio::EnsEMBL::Funcgen::Parsers::DataDumper;
my $parser = Bio::EnsEMBL::Funcgen::Parsers::DataDumper->new;

# This is what the parser will return. If it is not "used", the inheritance 
# in Bio::EnsEMBL::Funcgen::Probe to Bio::EnsEMBL::Storable will not work.
#
use Bio::EnsEMBL::Funcgen::Probe;

$parser->parse({
  data_dumper_file => $probe_features,
  call_back        => $map_transcript_to_genome,
});

$output_fh->close;

=head2 project_hit_to_genomic_coordinates

A hit looks like this:

$VAR1 = {
          'cigar_line' => '60=',
          'mismatch_count' => '0',
          'perc_id' => '100.00',
          'probe_seq_id' => '542',
          'q_end' => 60,
          'q_length' => 60,
          'q_start' => 0,
          'q_strand' => 1,
          'score' => '300',
          'scores' => 'scores:0:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:5:0:',
          't_end' => '1160',
          't_id' => 'ENST00000631077',
          't_length' => '1501',
          't_start' => 1101,
          't_strand' => 1,
          'total_mismatches' => 0
        };

  Transferred and adapted the code from the original ProbeAlign module. This
  could do with some refactoring.

=cut
sub project_hit_to_genomic_coordinates {

  my $param = shift;

  my $hit            = $param->{hit};
  my $genomic_blocks = $param->{genomic_blocks};
  my $transcript     = $param->{transcript};
  
  use Hash::Util qw( lock_hash );
  lock_hash(%$hit);
  
  my @trans_cigar_line = split ':', $hit->{cigar_line};

  my @stranded_cigar_line;
  if($transcript->strand == -1) {
    @$genomic_blocks      = reverse(@$genomic_blocks);
    @stranded_cigar_line  = reverse(@trans_cigar_line);
  }
  else{
    @stranded_cigar_line = @trans_cigar_line;
  }
  
  if(@$genomic_blocks==1) {
    my %projected_hit = %$hit;
    
    my $genomic_block = $genomic_blocks->[0];
    
    if ($genomic_block->isa('Bio::EnsEMBL::Mapper::Gap')) {
        warn(
            "This hit:\n"
            . Dumper($hit) . "\n"
            . "is projected into a gap:\n"
            . Dumper($genomic_block) . "\n"
        );
        return undef;
    }
    
    $projected_hit{t_start}  = $genomic_block->start;
    $projected_hit{t_end}    = $genomic_block->end;
    $projected_hit{t_strand} = $genomic_block->strand;
    
    $projected_hit{cigar_line} = join ' foo ', @stranded_cigar_line;
    
    return \%projected_hit;
  }

  #---------------------------------------------------------------------------
  #
  # Loop1
  #
  # This loop sets
  #
  # - $genomic_start and
  # - $genomic_end
  #
  # to the start and end positions of the probe feature in the genome.
  #
  # $gap_lengths{5} and $gap_lengths{3} are set to the number of bases
  # that the probe feature extends beyond the boundaries of the transcript
  # when mapped to the genome.
  #
  # @gaps is set to an array of pairs of start coordinates and lengths of 
  # gaps between transcript mapping locations.
  #
  my %gap_lengths = (
    5 => undef,
    3 => undef,
  );
  my @gaps;
  my $genomic_start;
  my $genomic_end;
  my $genomic_strand;
  
  GENOMIC_BLOCK: foreach my $block (@$genomic_blocks) {

    if(! $genomic_start) {
      if($block->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
      
        # All coordinate blocks will be on the same strand.
        $genomic_strand = $block->strand;
      
        # Set genomic_start
        if($gap_lengths{5}){
          # We have seen a gap
          $genomic_start = $block->start - $gap_lengths{5};
        } else {
          # Just set to first start
          $genomic_start = $block->start;
        }
      } else {
        #Must be 5' Gap/overhang
        $gap_lengths{5} = $block->length;
      }
    } elsif ($block->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
      $genomic_end = $block->end;
    } else {
      # Must be 3' Gap/overhang
      $gap_lengths{3}  = $block->length;
      $genomic_end    += $gap_lengths{3};
    }

    # Skip and Gap objects
    # as these will have already been inserted into the cigarline
    # These are 5'/3' overhangs
    next GENOMIC_BLOCK if $block->isa('Bio::EnsEMBL::Mapper::Gap');

    if (@gaps) {
      push @gaps, ($block->start - $gaps[$#gaps]);
    }
    push @gaps, ($block->end + 1);
  }

  #remove last value as this is the end of the match and not a gap start
  pop @gaps;
  #
  # End of loop1
  #---------------------------------------------------------------------------

  #---------------------------------------------------------------------------
  #
  # Loop2
  #
  # This loop deals with the cigar line at both the 5' and the 3' end.
  #

  # This is set to 5 and to 3.
  foreach my $end (keys %gap_lengths) {

    # If a part of the probe feature extends beyond the mapped transcript 
    # coordinates
    #
    if ($gap_lengths{$end}) {
    
      my $gap_block;

      if($end == 5) {
        $gap_block = shift @stranded_cigar_line;
      } else {
        $gap_block = pop @stranded_cigar_line;
      }

      # This parses the alignment block from the cigar line.
      my @tmp = split //, $gap_block;
      my $align_type = pop @tmp;
      (my $align_length = $gap_block) =~ s/$align_type//;

      if( $align_type ne 'X' ) {
        die("${end}' overhanging Gap has non-X/unexpected alignment type:\t$align_type");
      }

      my $new_gap_block = $gap_lengths{$end}.'S';
      
      # This seems to be trimming the match to the boundary of the mapped 
      # transcript.
      $align_length -= $gap_lengths{$end};

      if($align_length) {
        if($end == 5 ) {
          unshift @stranded_cigar_line, $align_length.$align_type;
        } else { 
          # is 3
          push @stranded_cigar_line, $align_length.$align_type;
        }
      }

      if($end == 5) {
        unshift @stranded_cigar_line, $new_gap_block;
      } else {
        push @stranded_cigar_line, $new_gap_block;
      }
    }
  }
  #
  # End of loop2
  #---------------------------------------------------------------------------

  #---------------------------------------------------------------------------
  #
  # Loop3
  #
  # This loop deals with the cigar line at both the 5' and the 3' end.
  #
  # Computes values for
  #
  my $match_count    = 0;
  my $mismatch_count = 0;
  my $score          = 0;
  my $q_length       = 0;
  my $cigar_line     = '';

  #Insert intron gaps into cDNA cigarline
  my $gap_start  = shift @gaps;
  my $gap_length = shift @gaps;

  #Set this to a hypothetical previous block end
  #So the initial block_start calculation will be valid
  my $block_end  = $genomic_start - 1;
  my $last_gap       = 0;
  
  foreach my $block (@stranded_cigar_line) {

    my @tmp = split //, $block;
    my $align_type = pop @tmp;
    (my $align_length = $block) =~ s/$align_type//;

    $q_length += $align_length;

    if ($align_type eq '=') {
      $match_count += $align_length;
      $score       += $align_length * 5;
    } else {
      # 'X' or S  mismatch (was U)
      $mismatch_count += $align_length;
      $score          -= $align_length * 4;
    }

    # We need to calculate the genomic end of the current block
    # given the previous block end or the genomic_start
    #
    my $block_start = $block_end + 1;
    $block_end += $align_length;

    if($block_end >= $gap_start) {

      #Could have multiple deletions
      while(($block_end >= $gap_start) && ! $last_gap) {
        #Insert the match first
        $align_length = $gap_start - $block_start;
        $cigar_line  .= $align_length.$align_type if $align_length;

        #Deletion wrt probe i.e. intron
        $cigar_line .= $gap_length.'D';

        #Now redefine start and end values
        #warn "block_start += $align_length + $gap_length";
        $block_start += $align_length + $gap_length;
        $block_end   += $gap_length;

        #Now grab the next gap
        if(@gaps) {
          $gap_start  = shift @gaps;
          $gap_length = shift @gaps;
        } else {
          $last_gap = 1;
        }
      }

      #We have reached the end of the gaps in this block
      #so just redefine block here
      $align_length = $block_end - $block_start + 1;
      $block = $align_length ? $align_length.$align_type : '';
    }
    $cigar_line .= $block;
  }
  #
  # End of loop3
  #---------------------------------------------------------------------------

  my %projected_hit = %$hit;

  $projected_hit{t_start}    = $genomic_start;
  $projected_hit{t_end}      = $genomic_end;
#   $projected_hit{strand}     = $transcript->strand;
  $projected_hit{t_strand}     = $genomic_strand;
  
  $projected_hit{cigar_line} = $cigar_line;
  
  return \%projected_hit;

#   # Test if we have already seen this alignment
#   my $gene = $gene_adaptor->fetch_by_transcript_stable_id($seq_id);
#   # The only way of doing this is to test the genomic_start/end and the genomic cigarline with the gene_stable_id and the probe_seq_id
#   $gene_sid = $gene->stable_id;
#   $gene_hit_key = "${gene_sid}:${probe_seq_id}:${genomic_start}:${genomic_end}:${cigar_line}";
# 
#   if(exists $gene_hits{$gene_hit_key}) {
#     $load_feature = 0;
#   } else {
#     #No need to count hits here
#     $gene_hits{$gene_hit_key} = undef;
#   }

#   # Now store the IDXref for this probe transcript hit
#   # This will mean we don't have recalculate this during the probe2transcript mapping step
#   #
#   $query_perc   = ($match_count/($match_count + $mismatch_count)) * 100;
#   $display_name = $self->get_display_name_by_stable_id($seq_id, 'transcript');
# 
# 
#   $xref = Bio::EnsEMBL::DBEntry->new (
#     -ANALYSIS => $feature->analysis,
#     -PRIMARY_ID => $seq_id,
#     -DISPLAY_ID => $display_name,
#     -DBNAME  => $edb_name,
#     -release => $schema_build,
#     -info_type => 'MISC',
#     -info_text => 'TRANSCRIPT',
#     -linkage_annotation => "ProbeTranscriptAlign $query_perc",#Add query_perc here when we have analysis
#     -version => $transcript_cache{$seq_id}->version, #version of transcript sid?
#   );
}

