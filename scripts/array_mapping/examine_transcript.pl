#!/usr/bin/env perl

use strict;
use locale; # Needed to ensure consistent ordering between Perl and UNIX sort
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Getopt::Long;
use Hash::Util qw( lock_hash lock_keys );

=head1

examine_transcript.pl \
  --registry /homes/mnuhn/work_dir_probemapping/lib/ensembl-funcgen/registry.pm \
  --species homo_sapiens \
  --transcript_probe_features_overlaps /nfs/nobackup/ensembl/mnuhn/probe2transcript/transcript_probe_features_overlaps.sorted_parallel.bed \
  --flanks_file         ./flanks.pl \
  --transcript_utr_file ./transcript_utr.pl \
  --transcript_info_file /nfs/nobackup/ensembl/mnuhn/probe2transcript/transcript_info.pl

=cut

my $registry;
my $species;
my $transcript_file;
my $flanks_file;
my $transcript_utr_file;
my $transcript_probe_features_overlaps;
my $transcript_info_file;
my $probe_feature_transcript_assignment_file;
my $probe_feature_transcript_rejection_file;

# Constants
my $debug = undef;
my $max_mismatches = 1;

GetOptions (
   'registry=s'                                  => \$registry,
   'species=s'                                   => \$species,
   'transcript_file=s'                           => \$transcript_file,
   'flanks_file=s'                               => \$flanks_file,
   'transcript_utr_file=s'                       => \$transcript_utr_file,
   'transcript_probe_features_overlaps=s'        => \$transcript_probe_features_overlaps,
   'transcript_info_file=s'                      => \$transcript_info_file,
   'probe_feature_transcript_assignments_file=s' => \$probe_feature_transcript_assignment_file,
   'probe_feature_transcript_rejection_file=s'   => \$probe_feature_transcript_rejection_file,
);

use Bio::EnsEMBL::Utils::Logger;
my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

Bio::EnsEMBL::Registry->load_all($registry);

use Bio::EnsEMBL::Funcgen::Parsers::DataDumper;
my $transcript_utr = Bio::EnsEMBL::Funcgen::Parsers::DataDumper->new->load_first_item_from_data_dump_file($transcript_utr_file);
my $flanks         = Bio::EnsEMBL::Funcgen::Parsers::DataDumper->new->load_first_item_from_data_dump_file($flanks_file);

use Hash::Util qw( lock_hash );
lock_hash(%$flanks);
$logger->info("We have " . (scalar keys %$flanks) . " flanks.\n");

my $funcgen_db_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
my $transcript_adaptor = Bio::EnsEMBL::Registry->get_adaptor  ($species, 'core', 'transcript');
my $all_core_transcripts = $transcript_adaptor->fetch_all();

open my $transcript_probe_features_overlaps_fh, '<', $transcript_probe_features_overlaps;
open my $transcript_info_fh , '>', $transcript_info_file;
open my $probe_feature_transcript_rejection_fh , '>', $probe_feature_transcript_rejection_file;

$logger->info("Writing to $probe_feature_transcript_assignment_file");

open my $probe_feature_transcript_assignment_fh , '>', $probe_feature_transcript_assignment_file;

# TODO turn HACK into something proper
sub add_xref {

  my $transcript_stable_id = shift;
  my $probe_feature_id     = shift;
  my $linkage_annotation   = shift;

  
  if (ref $linkage_annotation eq 'HASH') {
    use Carp;
    confess($linkage_annotation);
  }
  
  $probe_feature_transcript_assignment_fh->print(
    join "\t",
      $probe_feature_id,
      $transcript_stable_id,
      $linkage_annotation,
  );
  $probe_feature_transcript_assignment_fh->print("\n");
  return;
}

sub cache_and_load_unmapped_objects {

  my $param = shift;
  lock_hash(%$param);

  my $identifier       = $param->{identifier};
  my $object_type      = $param->{object_type};
  my $object_id        = $param->{object_id};
  my $summary          = $param->{summary};
  my $description      = $param->{description};
  
  my $unmapped_object_description = {
    dbID             => $object_id,
    object_type      => $object_type,
    stable_id        => $identifier,
    summary          => $summary,
    full_description => $description,
    type             => 'probe2transcript',
  };
  
  $probe_feature_transcript_rejection_fh->print(Dumper($unmapped_object_description));
  return;
}

associate_probes_to_transcripts({
  all_core_transcripts    => $all_core_transcripts,
  max_mismatches          => $max_mismatches,
  transcript_probe_features_overlaps_fh    => $transcript_probe_features_overlaps_fh,
  xref_db                 => $funcgen_db_adaptor,
  flanks                  => $flanks,
  transcript_info_file_fh => $transcript_info_fh,
});

$transcript_probe_features_overlaps_fh->close;
$transcript_info_fh->close;

$logger->finish_log;

=head2

  Description: OK this is where the action happens mostly. Goes through each overlap
  between transcript and probe feature, and collects data.
  Arg1: Array ref of Bio::EnsEMBL::Transcript objects
  Arg2: filehandle to the dump file created by bedtools
  Arg3: unmapped_counts: hashref
  Arg4: unmapped_objects: arrayref
  Arg5: xrefs: arrayref
  Arg6: Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor object
  Arg7: hashref with constant parameters
  Arg8: filehandle to log file
  Returntype: composite hashref: object_id -> transcript stable id -> arrayref -> 2 integers

=cut

sub associate_probes_to_transcripts {

  my $param = shift;
  
  my $all_core_transcripts    = $param->{all_core_transcripts};
  my $transcript_probe_features_overlaps_fh = $param->{transcript_probe_features_overlaps_fh};
  my $xref_db                 = $param->{xref_db};
  my $flanks                  = $param->{flanks};
  my $max_mismatches          = $param->{max_mismatches};
  my $transcript_info_file_fh = $param->{transcript_info_file_fh};
  
  my %all_core_transcripts_hashed_by_stable_id;
  foreach my $current_core_transcript (@$all_core_transcripts) {
  
    my $stable_id = $current_core_transcript->stable_id;
  
    die if (! defined $stable_id);
    die if (exists $all_core_transcripts_hashed_by_stable_id{$stable_id});
  
    $all_core_transcripts_hashed_by_stable_id{$current_core_transcript->stable_id} = $current_core_transcript;
  }

  examine_transcript({
    all_core_transcripts_hashed_by_stable_id => \%all_core_transcripts_hashed_by_stable_id,
    transcript_info_file_fh => $transcript_info_file_fh,
    max_mismatches          => $max_mismatches,
    pf_transc_overlap_fh    => $transcript_probe_features_overlaps_fh,
    xref_db                 => $xref_db,
    flanks                  => $flanks,
  });
}

=head2

  Description: For one trasncript, goes through each overlap
  between that transcript and probe feature, and collects data.
  CalledBy : associate_probes_to_transcripts
  Arg1: Bio::EnsEMBL::Transcript object
  Arg2: hashref to be filled
  Arg3: filehandle to the dump file created by bedtools
  Arg4: last line of the above file to be read
  Arg5: hashref for unmapped
  Arg6: arrayref of unmapped
  Arg7: arrayref of xrefs => create_final_xrefs
  Arg8: Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor object
  Arg9: hashref with constant parameters
  Arg10: filehandle to log file
  Returntype: String, last line read and first not processed (generally because it covers another transcript)

=cut

sub examine_transcript {

  my $param = shift;
  lock_hash(%$param);

  my $all_core_transcripts_hashed_by_stable_id = $param->{all_core_transcripts_hashed_by_stable_id};
  my $transcript_info_file_fh                  = $param->{transcript_info_file_fh};
  
  my $max_mismatches          = $param->{max_mismatches};
  my $pf_transc_overlap_fh    = $param->{pf_transc_overlap_fh};
  my $xref_db                 = $param->{xref_db};
  my $flanks                  = $param->{flanks};

  my $range_registry = Bio::EnsEMBL::Mapper::RangeRegistry->new();
  my $current_stable_id;
  
  my $unmapped_counts;
  my $unmapped_objects;
  my $xrefs;

  my $transcript_feature_info = {};
  
  # Process all following lines until you hit another transcript stable id
  OVERLAP_LINE: while (defined (my $line = <$pf_transc_overlap_fh>)) {

    chomp $line;
    
    if ($debug) {
      $logger->info("LINE:$line\n");
    }
    
    my (
      $tx_chr, 
      $tx_start, 
      $tx_end, 
      $transcript_sid_2, 
      $chr, 
      $start, 
      $end, 
      $strand, 
      $cigar, 
      $mismatches, 
      $feature_id, 
      $probe_id, 
      $probe_name, 
      $probeset_id, 
      $probeset_name, 
      $array_name, 
      $array_vendor, 
      $array_class,
      $is_probeset_array,
      $is_linked_array,
      $has_sense_interrogation
    ) = split "\t", $line;
    #     1       11869    14409     ENST00000456328    1     11869   11893   1       25=     0           22771302     10249677     578754    1631966        8044649 HuGene-1_0-st-v1        AFFY
    
    die($transcript_sid_2) unless(exists $all_core_transcripts_hashed_by_stable_id->{$transcript_sid_2});
    my $transcript = $all_core_transcripts_hashed_by_stable_id->{$transcript_sid_2};
    
    # Check, if the current line is the start of overlaps from the next 
    # transcript.
    #
    if ($current_stable_id ne $transcript_sid_2) {

      # $current_stable_id will be undefined in the first iteration. Only
      # print results in the following iterations.
      #
      if (defined $current_stable_id) {
        $transcript_info_file_fh->print(
          Dumper({
            transcript_stable_id => $current_stable_id,
            probe_hits_by_array  => $transcript_feature_info,
          })
        );
        $transcript_feature_info = {};
      }
    
      $current_stable_id = $transcript_sid_2;
    
      my @exons = @{$transcript->get_all_Exons};
      $range_registry->flush();
      
      foreach my $exon (@exons) {
        $range_registry->check_and_register('exonic', $exon->seq_region_start, $exon->seq_region_end);
      }
      my $first_exon = $exons[0];
      my $last_exon  = $exons[$#exons];

      my %exonutrs;
      my $transcript_sid = $transcript->stable_id;

      $logger->debug("TRANSCRIPT $transcript_sid\n");

      my $transcript_slice = $transcript->feature_Slice();
      my $slice            = $transcript_slice->expand($flanks->{$transcript_sid}{5}, $flanks->{$transcript_sid}{3});
      # Find flanking regions
      if ($transcript->strand == 1) {
        if ($flanks->{$transcript_sid}{3}) {
          $exonutrs{3} = [$last_exon->seq_region_start, $slice->end];
        }
        if ($flanks->{$transcript_sid}{5}) {
          $exonutrs{5} = [$slice->start, $first_exon->seq_region_end];
        }
      } else {
        if ($flanks->{$transcript_sid}{3}) {
          $exonutrs{3} = [$slice->start(), $last_exon->seq_region_end];
        }
        if ($flanks->{$transcript_sid}{5}) {
          $exonutrs{5} = [$first_exon->seq_region_start, $slice->end];
        }
      }

      # Mark flanking regions
      foreach my $end(3, 5) {
        if ($flanks->{$transcript_sid}{$end}) {
          $range_registry->check_and_register("${end}_exonutr", @{$exonutrs{$end}});
        }
      }
    }

    examine_probefeature({
      transcript              => $transcript,
      feature_id              => $feature_id,
      probe_id                => $probe_id,
      probe_name              => $probe_name,
      probeset_id             => $probeset_id,
      probeset_name           => $probeset_name,
      start                   => $start,
      end                     => $end,
      strand                  => $strand,
      cigar                   => $cigar,
      mismatches              => $mismatches,
      transcript_feature_info => $transcript_feature_info,
      unmapped_counts         => $unmapped_counts,
      unmapped_objects        => $unmapped_objects,
      xrefs                   => $xrefs,
      xref_db                 => $xref_db,
      range_registry          => $range_registry,
      array_name              => $array_name, 
      array_vendor            => $array_vendor,
      array_class             => $array_class,

      is_probeset_array       => $is_probeset_array,
      has_sense_interrogation => $has_sense_interrogation,

      flanks                  => $flanks,
    });
  }
  $transcript_info_file_fh->print(
    Dumper({
      transcript_stable_id => $current_stable_id,
      probe_hits_by_array  => $transcript_feature_info,
    })
  );

  return;
}

=head2

  Description: For one overlap between a trasncript and a probe feature, collects data
  CalledBy: examine_transcript()
  Arg1: Bio::EnsEMBL::Transcript object
  Arg2: ProbeFeature DB id
  Arg3: Probe DB id
  Arg4: ProbeSet DB id
  Arg5: ProbeSet name, string
  Arg6: start of probe
  Arg7: end of probe
  Arg8: strand (1 or -1)
  Arg9: CIGAR string of probefeature alignment
  Arg10: number of mismatches
  Arg11: hashref to be filled
  Arg12: hashref of stats on unmapped objects
  Arg13: arrayref of Bio::EnsEMBL::UnmappedObject objects
  Arg14: arrayref for add_xref()
  Arg15: Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor object
  Arg16: hashref with constant parameters
  Arg17: filehandle to log file

=cut

sub examine_probefeature {

  my $param = shift;  
  lock_hash(%$param);

  my $transcript              = $param->{transcript};
  my $feature_id              = $param->{feature_id};
  my $probe_id                = $param->{probe_id};
  my $probe_name              = $param->{probe_name};
  my $start                   = $param->{start};
  my $end                     = $param->{end};
  my $strand                  = $param->{strand};
  my $cigar_line              = $param->{cigar};
  my $mismatch_count          = $param->{mismatches};
  my $transcript_feature_info = $param->{transcript_feature_info};
  my $unmapped_counts         = $param->{unmapped_counts};
  my $unmapped_objects        = $param->{unmapped_objects};
  my $xrefs                   = $param->{xrefs};
  my $xref_db                 = $param->{xref_db};
  my $array_name              = $param->{array_name};
  my $array_vendor            = $param->{array_vendor};
  my $array_class             = $param->{array_class};
  my $probeset_name           = $param->{probeset_name};
  my $probeset_id             = $param->{probeset_id};
  my $range_registry          = $param->{range_registry};
  my $flanks                  = $param->{flanks};

  my $is_probeset_array       = $param->{is_probeset_array};
  my $has_sense_interrogation = $param->{has_sense_interrogation};

  my $transcript_sid = $transcript->stable_id();
  my $transcript_version     = $transcript->version;
  my $log_name;

  if ($is_probeset_array) {
    $log_name       = $transcript_sid."\t(${probeset_name})\t${probe_id}";
  } else{
    $log_name       = $transcript_sid."\t(".$probe_name.")\t${probe_id}";
  }

#   The ultimate guide to sense interrogation
#   =========================================
# 
#   This is how sense interrogation is meant to play out:
# 
#             Does the array have sense interrogation?
# 
#                               |
#                               |
#                 ----------------------------
#                 |                          |
#                 |                          |
# 
#                 No                         Yes
# 
#                 |                          |
#                 |                          |
# 
#           Is the probe feature        Is the probe feature 
#           on the same strand          on the same strand  
#           as the transcript?          as the transcript?
# 
#                 |                           |
#                 |                           |
#           --------------            ---------------
#           |            |            |             |
#           |            |            |             |
#
#           Yes          No           Yes           No
#
#           |            |            |             |
#           |            |            |             |
#         Accept        Reject       Reject        Accept   
#         the match     the match    the match     the match
#
#   Note that for probes from a sense interrogation array probe features 
#   are accepted only, if the sequence matches on the reverse strand of 
#   the transcript.

  if ($has_sense_interrogation) {
    if($transcript->seq_region_strand == $strand) {
      if ($debug) {
        $logger->info('Unmapped sense '.$log_name."\n");
      }
      cache_and_load_unmapped_objects({
        identifier  => $transcript_sid,
        object_type => 'ProbeFeature',
        object_id   => $feature_id,
        summary     => 'Unmapped sense',
        description => 'Unmapped sense '.$log_name,
      });
      return;
    }
  } elsif ($transcript->seq_region_strand != $strand) {
    if ($debug) {
      $logger->info('Unmapped anti-sense '.$log_name."\n");
    }
    cache_and_load_unmapped_objects({
      identifier  => $transcript_sid,
      object_type => 'ProbeFeature',
      object_id   => $feature_id,
      summary     => 'Unmapped anti-sense',
      description => 'Unmapped anti-sense '.$log_name,
    });

    return;
  }

  my $mm_link_txt = '';
  my $has_mismatches = 0;
  if ($mismatch_count) {
    $has_mismatches = 1;
    $mm_link_txt = ' ('. $mismatch_count.' bp mismatch)';
  }

  if($cigar_line =~ /D/) {
    record_gapped_probefeature({
      feature_id              => $feature_id,
      probe_id                => $probe_id,
      probeset_id             => $probeset_id,
      transcript              => $transcript, 
      transcript_feature_info => $transcript_feature_info,
      mm_link_txt             => $mm_link_txt,
      has_mismatches          => $has_mismatches, 
      log_name                => $log_name,
      unmapped_counts         => $unmapped_counts,
      unmapped_objects        => $unmapped_objects,
      xref_db                 => $xref_db,
      array_name              => $array_name,

      is_probeset_array       => $is_probeset_array,
#       is_linked_array         => $is_linked_array,
      has_sense_interrogation => $has_sense_interrogation,

    });
  } else {
    record_aligned_probefeature({
      feature_id              => $feature_id,
      probe_id                => $probe_id, 
      probeset_id             => $probeset_id, 
      start                   => $start, 
      end                     => $end, 
      cigar_line              => $cigar_line, 
      mismatch_count          => $mismatch_count,
      transcript_feature_info => $transcript_feature_info,
      range_registry          => $range_registry,
      transcript              => $transcript,
      mm_link_txt             => $mm_link_txt,
      has_mismatches          => $has_mismatches, 
      log_name                => $log_name, 
      unmapped_counts         => $unmapped_counts, 
      unmapped_objects        => $unmapped_objects, 
      xrefs                   => $xrefs, 
      xref_db                 => $xref_db, 
      array_name              => $array_name,

      is_probeset_array       => $is_probeset_array,
      has_sense_interrogation => $has_sense_interrogation,

      flanks                  => $flanks,
    });
  }
  return;
}

=head2

  Description: For one gapped overlap between a trasncript and a probe feature, collects data
  
  The data is returned by making changes to the hash referenced by 
  $transcript_feature_info instead of a return value.
  
  CalledBy: examine_probefeature
  Arg1: ProbeFeature DB id
  Arg2: Probe DB id
  Arg3: ProbeSet DB id
  Arg4: Bio::EnsEMBL::Transcript object
  Arg5: hashref to fill
  Arg6: Text string describing overlap
  Arg7: boolean (contains mismatches?)
  Arg8: string describing the overlapped pair
  Arg9: hashref of stats on unmapped objects
  Arg10: arrayref of Bio::EnsEMBL::UnmappedObject objects
  Arg11: hashref with constant parameters
  Arg12: filehandle to log file

=cut

sub record_gapped_probefeature {
  
  my $param = shift;
  lock_hash(%$param);
  
  my $feature_id              = $param->{feature_id};
  my $probe_id                = $param->{probe_id};
  my $probeset_id             = $param->{probeset_id};
  my $transcript              = $param->{transcript};
  my $transcript_feature_info = $param->{transcript_feature_info};
  my $mm_link_txt             = $param->{mm_link_txt};
  my $has_mismatches          = $param->{has_mismatches};
  my $log_name                = $param->{log_name};
  my $unmapped_counts         = $param->{unmapped_counts};
  my $unmapped_objects        = $param->{unmapped_objects};
  my $xref_db                 = $param->{xref_db};
  my $array_name              = $param->{array_name};

  my $is_probeset_array       = $param->{is_probeset_array};
#   my $is_linked_array         = $param->{is_linked_array};
  my $has_sense_interrogation = $param->{has_sense_interrogation};

  my $transcript_sid = $transcript->stable_id();

  my $sql = 'select count(hit_id) as num_transcripts from probe_feature where probe_feature_id=? and hit_id=? and source="transcript"';
  
  my $sql_helper = $xref_db->dbc->sql_helper;
  my $num_transcripts = $sql_helper->execute_single_result(
    -SQL    => $sql,
    -PARAMS => [ 
      $feature_id, 
      $transcript_sid 
    ]
  );
  
  if($num_transcripts==0) {
    if ($debug) {
      $logger->info('Unmapped Gapped ProbeFeature '.$log_name."\n");
    }
    cache_and_load_unmapped_objects({
      identifier  => $transcript_sid,
      object_type => 'ProbeFeature',
      object_id   => $feature_id,
      summary     => 'Unmapped Gapped ProbeFeature',
      description => 'Gapped ProbeFeature did not match transcript structure',
    });
    return
  }

  if ($debug) {
    $logger->info('Mapped Gapped ProbeFeature '.$log_name."\n");
  }
  my $probe_match_annotation = {
    annotation => "exon-exon match${mm_link_txt}", 
    has_mismatches => $has_mismatches,
  };

  if ($is_probeset_array) {
  
    if (! exists $transcript_feature_info->{$array_name}) {
    
      # The structure for probesets:
      $transcript_feature_info->{$array_name} = {
        probeset_array => 1,
        array_name     => $array_name,
        probesets      => {},
      };
      lock_keys(%{$transcript_feature_info->{$array_name}});
    }

    if (! exists $transcript_feature_info->{$array_name}->{probesets}->{$probeset_id}) {
      $transcript_feature_info->{$array_name}->{probesets}->{$probeset_id} = {
        probe_id                   => {},
        num_probes_with_mismatches => 0,
        num_probe_features         => 0,
      };
      lock_keys(%{$transcript_feature_info->{$array_name}->{probesets}->{$probeset_id}});
    }
    if (! exists $transcript_feature_info->{$array_name}->{probesets}->{$probeset_id}->{probe_id}->{$probe_id}) {
      $transcript_feature_info->{$array_name}->{probesets}->{$probeset_id}->{probe_id}->{$probe_id} = [];
    }
    
    $transcript_feature_info->{$array_name}->{probesets}->{$probeset_id}->{num_probe_features} += 1;
    if ($has_mismatches) {
      $transcript_feature_info->{$array_name}->{probesets}->{$probeset_id}->{num_probes_with_mismatches} += 1;
    }
    push 
      @{$transcript_feature_info->{$array_name}->{probesets}->{$probeset_id}->{probe_id}->{$probe_id}}, 
      $probe_match_annotation;

  } else {

    if (! exists $transcript_feature_info->{$array_name}) {
      #
      # The structure for arrays without probesets:
      #
      $transcript_feature_info->{$array_name} = {
        probeset_array => undef,
        array_name     => $array_name,
        probe          => {},
      };
      lock_keys(%{$transcript_feature_info->{$array_name}});
    }
    if (! exists $transcript_feature_info->{$array_name}->{probe}->{$probe_id}) {
      $transcript_feature_info->{$array_name}->{probe}->{$probe_id} = [];
    }
    push 
      @{$transcript_feature_info->{$array_name}->{probe}->{$probe_id}}, 
      $probe_match_annotation;
  }
  add_xref($transcript_sid, $feature_id, $probe_match_annotation->{annotation});
  # The real return value is in $transcript_feature_info and data is written in cache_and_load_unmapped_objects
  return;
}

=head2

  Description: For one ungapped overlap between a trasncript and a probe feature, collects data
  
  The data is returned by making changes to the hash referenced by 
  $transcript_feature_info instead of a return value.
  
  It also stores xrefs.
  
  CalledBy: examine_probefeature()
  Arg1: ProbeFeature DB id
  Arg2: Probe DB id
  Arg3: ProbeSet DB id
  Arg4: start of probe
  Arg5: end of probe
  Arg6: CIGAR string
  Arg7: integer count of mismatches
  Arg8: hashref to fill
  Arg9: Bio::EnsEMBL::Transcript object
  Arg10: Text string describing overlap
  Arg11: boolean if contains mismatches
  Arg12: string describing the overlapped pair
  Arg13: hashref of stats on unmapped objects
  Arg14: arrayref of Bio::EnsEMBL::UnmappedObject objects
  Arg15: arrayref for add_xref()
  Arg16: Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor object
  Arg17: hashref with constant parameters
  Arg18: filehandle to log file

=cut

sub record_aligned_probefeature {
  my $param = shift;
  lock_hash(%$param);

  my $feature_id              = $param->{feature_id};
  my $probe_id                = $param->{probe_id};
  my $probeset_id             = $param->{probeset_id};
  my $start                   = $param->{start};
  my $end                     = $param->{end};
  my $cigar_line              = $param->{cigar_line};
  my $mismatch_count          = $param->{mismatch_count};
  my $has_mismatches          = $param->{has_mismatches};
  my $transcript_feature_info = $param->{transcript_feature_info};
  my $transcript              = $param->{transcript};
  my $mm_link_txt             = $param->{mm_link_txt};
  my $log_name                = $param->{log_name};
  my $unmapped_counts         = $param->{unmapped_counts};
  my $unmapped_objects        = $param->{unmapped_objects};
  my $xrefs                   = $param->{xrefs};
  my $xref_db                 = $param->{xref_db};
  my $range_registry          = $param->{range_registry};
  my $flanks                  = $param->{flanks};

  my $array_name              = $param->{array_name};

  my $is_probeset_array       = $param->{is_probeset_array};
  my $has_sense_interrogation = $param->{has_sense_interrogation};

  my $five_mismatch  = 0;
  my $three_mismatch = 0;
  my $feature_start  = $start;
  my $feature_end    = $end;
  
  if($cigar_line =~ /(^[0-9]+)m/) {
    $five_mismatch = $1;
    $feature_start += $five_mismatch;
  }

  if($cigar_line =~ /([0-9]+$)m/) {
    $three_mismatch = $1;
    $feature_end -= $three_mismatch;
  }

  my $min_overlap  = ($end - $start + 1 - $max_mismatches + ($mismatch_count - $five_mismatch - $three_mismatch));
  my $exon_overlap = $range_registry->overlap_size('exonic', $feature_start, $feature_end);
  my $flank_end     = 0;
  my $flank_overlap = 0;

  foreach my $side (3, 5) {
    if ($flanks->{$transcript->stable_id}->{$side}) {
      $flank_overlap = $range_registry->overlap_size("${side}_exonutr", $start, $end);
    }
    if ($flank_overlap) {
      $flank_end = $side;
      last;
    }
  }
  
  my $overlap_is_good_enough = ($exon_overlap >= $min_overlap) || ($flank_overlap >= $min_overlap);

  if ($overlap_is_good_enough) {
    my $linkage_annotation;

    if ($exon_overlap && $flank_overlap) {
      $linkage_annotation = "exon/${flank_end}' flank boundary${mm_link_txt}";
    } elsif ($exon_overlap) {
      $linkage_annotation = "exon${mm_link_txt}";
    } else {                                #only flank over lap
      $linkage_annotation = "${flank_end}' flank${mm_link_txt}";
    }
    my $probe_match_annotation = {
      annotation     => $linkage_annotation, 
      has_mismatches => $has_mismatches,
    };

    if ($is_probeset_array) {
    
      # This is the output of the function
      if (! exists $transcript_feature_info->{$array_name}) {
      
        # The structure for probesets:
        $transcript_feature_info->{$array_name} = {
          probeset_array => 1,
          array_name     => $array_name,
          probesets      => {},
        };
        lock_keys(%{$transcript_feature_info->{$array_name}});
      }

      if (! exists $transcript_feature_info->{$array_name}->{probesets}->{$probeset_id}) {
        $transcript_feature_info->{$array_name}->{probesets}->{$probeset_id} = {
          probe_id           => {},
          num_mismatches     => 0,
          
          # Number of probe features from this probe set that match on this transcript.
          num_probe_features => 0,
        };
        lock_keys(%{$transcript_feature_info->{$array_name}->{probesets}->{$probeset_id}});
      }
      if (! exists $transcript_feature_info->{$array_name}->{probesets}->{$probeset_id}->{probe_id}->{$probe_id}) {
        $transcript_feature_info->{$array_name}->{probesets}->{$probeset_id}->{probe_id}->{$probe_id} = [];
      }
      
      $transcript_feature_info->{$array_name}->{probesets}->{$probeset_id}->{num_probe_features} += 1;
      if ($has_mismatches) {
        $transcript_feature_info->{$array_name}->{probesets}->{$probeset_id}->{num_mismatches} += 1;
      }
      push 
        @{$transcript_feature_info->{$array_name}->{probesets}->{$probeset_id}->{probe_id}->{$probe_id}},
        $probe_match_annotation;
      
    } else {
    
      if (! exists $transcript_feature_info->{$array_name}) {
        #
        # The structure for arrays without probesets:
        #
        $transcript_feature_info->{$array_name} = {
          probeset_array => undef,
          array_name     => $array_name,
          probe          => {},
        };
        lock_keys(%{$transcript_feature_info->{$array_name}});
      }
      if (! exists $transcript_feature_info->{$array_name}->{probe}->{$probe_id}) {
        $transcript_feature_info->{$array_name}->{probe}->{$probe_id} = [];
      }
      push 
        @{$transcript_feature_info->{$array_name}->{probe}->{$probe_id}}, 
        $probe_match_annotation;
    }
    if ($debug) {
      $logger->info("Mapped\t$probe_id\t$probeset_id\t".$transcript->stable_id()."\n");
    }
    add_xref($transcript->stable_id, $feature_id, $linkage_annotation);
  } else {
    my ($summary, $region);

    if (!( $exon_overlap || $flank_overlap)) {
      $summary = 'intronic';
      $region  = 'intronic region';
    } elsif ($exon_overlap) {
      if (! $flank_overlap) {
        $summary = 'exon boundary';
        $region  = $summary;
      } else {
        $summary = "${flank_end}' flank boundary";
        $region  = $summary;
      }
    } else {
      $summary = "${flank_end}' flank boundary";
      $region  = $summary;
    }
    if ($debug) {
      $logger->info("Unmapped $summary ".$log_name."\n");
    }
    cache_and_load_unmapped_objects({
      identifier  => $transcript->stable_id,
      object_type => 'ProbeFeature',
      object_id   => $feature_id,
      summary     => "Unmapped $summary",
      description => "Probe mapped to $region of transcript",
    });
  }
  # The real return value is in $transcript_feature_info
  return;
}

