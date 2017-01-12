#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Getopt::Long;
use Hash::Util qw( lock_hash );


=head1

bsub -M4000 -R"select[mem>4000] rusage[mem=4000]" -Is $SHELL

# Takes about 90 minutes
# 23 minutes upon rerun and 3.3GB ram
# 
calculate_utrs.pl \
  --registry /homes/mnuhn/work_dir_probemapping/lib/ensembl-funcgen/registry.pm \
  --species homo_sapiens \
  --transcript_utr_file ./unannotated_utrs.pl

# 32 minutes
# 4Gb Ram

write_extended_transcripts_into_file.pl \
  --registry /homes/mnuhn/work_dir_probemapping/lib/ensembl-funcgen/registry.pm \
  --species homo_sapiens \
  --extended_transcript_outputfile ./transcripts.bed \
  --unannotated_utrs ./unannotated_utrs.pl \
  --flanks_outputfile ./flanks.pl

bedSort ./transcripts.bed ./transcripts.bed

# Takes 23 minutes and generates a 45GB file:

date
mysql --quick -NB --host=mysql-ens-reg-prod-1 --port=4526 --user=ensro mnuhn_alnnew_homo_sapiens_funcgen_86_38 -e "
  select
    seq_region.name, 
    seq_region_start, 
    seq_region_end, 
    seq_region_strand, 
    cigar_line, 
    mismatches, 
    probe_feature_id, 
    probe_id, 
    probe.
    name, 
    probe_set_id, 
    probe_set.name,
    array.name,
    array.vendor,
    array.class
  from
    probe_feature join probe using(probe_id)
    left join probe_set using(probe_set_id)
    join array_chip using(array_chip_id)
    join array using(array_id)
    join seq_region using(seq_region_id)
" > /nfs/nobackup/ensembl/mnuhn/probe2transcript/probe_features_nogroup.bed
date

# Uses 6.7 GB ram:

# bedSort /nfs/nobackup/ensembl/mnuhn/probe2transcript/probe_features.bed /nfs/nobackup/ensembl/mnuhn/probe2transcript/probe_features.bed
# bedSort /nfs/nobackup/ensembl/mnuhn/probe2transcript/probe_features_nogroup.bed /nfs/nobackup/ensembl/mnuhn/probe2transcript/probe_features_nogroup.bed

wc -l /nfs/nobackup/ensembl/mnuhn/probe2transcript/probe_features_nogroup.bed
489291541 /nfs/nobackup/ensembl/mnuhn/probe2transcript/probe_features_nogroup.bed

echo "489291541 / 500" | bc
978583

mkdir -p /nfs/nobackup/ensembl/mnuhn/probe2transcript/probe_features_nogroup_split/

split --lines=978583 /nfs/nobackup/ensembl/mnuhn/probe2transcript/probe_features_nogroup.bed /nfs/nobackup/ensembl/mnuhn/probe2transcript/probe_features_nogroup_split/

split_files=`find /nfs/nobackup/ensembl/mnuhn/probe2transcript/probe_features_nogroup_split/ -type f`

for split_file in $split_files
do
  bsub bedSort $split_file ${split_file}.sorted
done

sorted_split_files=`find /nfs/nobackup/ensembl/mnuhn/probe2transcript/probe_features_nogroup_split/ -type f | grep sorted | sort`

date
sort --unique --buffer-size=2G -T /nfs/nobackup/ensembl/mnuhn/probe2transcript/ -m $sorted_split_files > /nfs/nobackup/ensembl/mnuhn/probe2transcript/probe_features_nogroup.sorted.bed
date

date
bedSort /nfs/nobackup/ensembl/mnuhn/probe2transcript/probe_features_nogroup_split/cs /nfs/nobackup/ensembl/mnuhn/probe2transcript/probe_features_nogroup_split/cs.sorted
date

# This takes 9 minutes

date
sort --unique -k 1,1n -k 2,2n -k 3,3n --buffer-size=2G -T /nfs/nobackup/ensembl/mnuhn/probe2transcript/ /nfs/nobackup/ensembl/mnuhn/probe2transcript/probe_features_nogroup.bed > /nfs/nobackup/ensembl/mnuhn/probe2transcript/probe_features_nogroup.sorted.bed
date

# sort sorts, but the numerical sort doesn't seem to quite work
# Needs 8GB ram
bedSort /nfs/nobackup/ensembl/mnuhn/probe2transcript/probe_features_nogroup.sorted.bed /nfs/nobackup/ensembl/mnuhn/probe2transcript/probe_features_nogroup.sorted.bed

#  5 minutes
date
bedtools intersect -sorted -wa -wb -a ./transcripts.bed -b /nfs/nobackup/ensembl/mnuhn/probe2transcript/probe_features_nogroup.sorted.bed > /nfs/nobackup/ensembl/mnuhn/probe2transcript/transcript_probe_features_overlaps.bed
date




bsub -n 16 -q production-rh7 -M32000 -R"select[mem>32000] rusage[mem=32000]" -Is $SHELL
date
sort --parallel=16 -u -T /nfs/nobackup/ensembl/mnuhn/probe2transcript/ --buffer-size=31G -k4,4 /nfs/nobackup/ensembl/mnuhn/probe2transcript/probe_features_nogroup.test.bed > /nfs/nobackup/ensembl/mnuhn/probe2transcript/probe_features_nogroup.test.sorted.bed
date

bsub -n 32 -q production-rh7 -M65000 -R"select[mem>65000] rusage[mem=65000]" -Is $SHELL
date
sort --parallel=16 -T /nfs/nobackup/ensembl/mnuhn/probe2transcript/ --buffer-size=64G -k4,4 /nfs/nobackup/ensembl/mnuhn/probe2transcript/probe_features_nogroup.test.bed > /nfs/nobackup/ensembl/mnuhn/probe2transcript/probe_features_nogroup.test.sorted.bed
date





bsub -n 16 -q production-rh7 -M32000 -R"select[mem>32000] rusage[mem=32000]" -Is $SHELL
date
sort --parallel=16 -T /nfs/nobackup/ensembl/mnuhn/probe2transcript/ --buffer-size=30G -k4,4 /nfs/nobackup/ensembl/mnuhn/probe2transcript/transcript_probe_features_overlaps.bed > /nfs/nobackup/ensembl/mnuhn/probe2transcript/transcript_probe_features_overlaps.sorted_parallel.bed
date

# # Slow!
# date
# sort -T /nfs/nobackup/ensembl/mnuhn/probe2transcript/ --buffer-size=512M -k4,4 /nfs/nobackup/ensembl/mnuhn/probe2transcript/transcript_probe_features_overlaps.bed > /nfs/nobackup/ensembl/mnuhn/probe2transcript/transcript_probe_features_overlaps.sorted.bed
# date

calculate_arrays_per_object.pl \
  --registry /homes/mnuhn/work_dir_probemapping/lib/ensembl-funcgen/registry.pm \
  --species homo_sapiens \
  --arrays_per_object_file /nfs/nobackup/ensembl/mnuhn/probe2transcript/arrays_per_object.pl \
  --probeset_sizes_file    /nfs/nobackup/ensembl/mnuhn/probe2transcript/probeset_sizes.pl    \
  --object_names_file      /nfs/nobackup/ensembl/mnuhn/probe2transcript/object_names.pl      \

examine_transcript.pl \
  --registry /homes/mnuhn/work_dir_probemapping/lib/ensembl-funcgen/registry.pm \
  --species homo_sapiens \
  --transcript_probe_features_overlaps /nfs/nobackup/ensembl/mnuhn/probe2transcript/transcript_probe_features_overlaps.sorted_parallel.bed \
  --flanks_file         ./flanks.pl \
  --transcript_utr_file ./transcript_utr.pl \
  --transcript_info_file /nfs/nobackup/ensembl/mnuhn/probe2transcript/transcript_info.pl

compute_hits.pl \
  --probeset_sizes_file   /nfs/nobackup/ensembl/mnuhn/probe2transcript/probeset_sizes.pl   \
  --transcript_info_file  /nfs/nobackup/ensembl/mnuhn/probe2transcript/transcript_info.pl \
  --probeset_transcript_assignments /nfs/nobackup/ensembl/mnuhn/probe2transcript/probeset_transcript_assignments.pl \
  --probeset_transcript_rejections  /nfs/nobackup/ensembl/mnuhn/probe2transcript/probeset_transcript_rejections.pl \
  --probe_transcript_assignments    /nfs/nobackup/ensembl/mnuhn/probe2transcript/probe_transcript_assignments.pl \


=cut

=head1

  Writes coordinates of extended transcripts to an unsorted bed file.

=cut

# Input
my $registry;
my $species;
my $unannotated_utrs_file;

# The names of the output files
my $flanks_outputfile;
my $extended_transcript_outputfile;

GetOptions (
   'registry=s'         => \$registry,
   'species=s'          => \$species,
   'unannotated_utrs=s' => \$unannotated_utrs_file,
   
   'extended_transcript_outputfile=s' => \$extended_transcript_outputfile,
   'flanks_outputfile=s'              => \$flanks_outputfile,
);

use Bio::EnsEMBL::Utils::Logger;
my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

use Bio::EnsEMBL::Funcgen::Parsers::DataDumper;
my $unannotated_utrs = Bio::EnsEMBL::Funcgen::Parsers::DataDumper->new->load_first_item_from_data_dump_file($unannotated_utrs_file);
lock_hash(%$unannotated_utrs);

# Not currently used, just transferred this from probe2transcript.
my $utr_multiplier;

# Also not currently used, just transferred this from probe2transcript.
my $utr_extends = {
  3 => undef,
  5 => undef,
};

assert_parameters_ok({
  utr_multiplier   => $utr_multiplier,
  utr_extends      => $utr_extends,
  unannotated_utrs => $unannotated_utrs,
});

Bio::EnsEMBL::Registry->load_all($registry);

$logger->info("Fetching all transcripts from core database\n");

my $transcript_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'transcript');
my $transcripts        = $transcript_adaptor->fetch_all;

$logger->info("Done.\n");

$logger->info("Starting write_extended_transcripts_into_file\n");

(
  my $flanks, 
  my $utr_counts
) = write_extended_transcripts_into_file({
  transcripts      => $transcripts, 
  unannotated_utrs => $unannotated_utrs,
  utr_multiplier   => $utr_multiplier,
  utr_extends      => $utr_extends,
  
  extended_transcript_outputfile  => $extended_transcript_outputfile,
});

$logger->info("Done write_extended_transcripts_into_file\n");

# Prints out some random stats
foreach my $end(5, 3) {
  my $total_length = 0;
  foreach my $utr_count (@{$utr_counts->{$end}}) {
    $total_length += $utr_count;
  }
  my $num_utrs = scalar(@{$utr_counts->{$end}});
  my $average = ($num_utrs) ? ($total_length/$num_utrs) : 0;
  $logger->info("Seen $num_utrs $end prime UTRs with and average length of $average\n");
}

$logger->info("Writing results to $flanks_outputfile\n");

open my $flanks_outputfile_fh, '>', $flanks_outputfile;
$flanks_outputfile_fh->print(Dumper($flanks));
$flanks_outputfile_fh->close;

$logger->info("Done writing results.\n");

$logger->finish_log;

sub assert_parameters_ok {
  my $param = shift;
  
  my $utr_multiplier   = $param->{utr_multiplier};
  my $utr_extends      = $param->{utr_extends};
  my $unannotated_utrs = $param->{unannotated_utrs};
  
  #Multiplier needs to be a pisitive real number
  if(defined $utr_multiplier && $utr_multiplier !~ /-*[0-9]+[\.0-9]*/) {
    $logger->error("-utr_multiplier must be a positive real number:\t$utr_multiplier" . "\n");
  }

  #Validate extend params
  for my $end (3, 5) {
    if(defined $utr_extends->{$end} && $utr_extends->{$end} =~ /\D+/) {
      $logger->error("Invalid -${end}_prime_extend parameter(".$utr_extends->{$end}.").  Must be a number(bp)" . "\n");
    }
  }

  #Validate unannotated defaults
  for my $end (3, 5) {
    if(defined $unannotated_utrs->{$end}) {
      if($unannotated_utrs->{$end} =~ /^\D+$/) {
        $logger->error("Invalid -unannotated_${end}_utr parameter(".$unannotated_utrs->{$end}.").  Must be a number(bp)" . "\n");
      }
      else{
        $logger->info("Setting ${end} unannotated UTR length to ".$unannotated_utrs->{$end} . "\n");
      }
    }
    else{
      if(! $unannotated_utrs) {
        $logger->info("Defaulting unannotated ${end} UTR length to 0" . "\n");
        $unannotated_utrs->{$end} = 0;
      }
    }
  }
}

=head2 write_extended_transcripts_into_file

  Description: Dumps transcript regions into sorted BED file

  Arg1: Array ref of Bio::EnsEMBL::Transcript objects
  Arg2: Filename to be written into
  Arg3: Hash ref to collect stats

=cut

sub write_extended_transcripts_into_file {

  my $param = shift;
  lock_hash(%$param);

  my $transcripts      = $param->{transcripts};
  my $filename         = $param->{extended_transcript_outputfile};
  my $unannotated_utrs = $param->{unannotated_utrs};
  my $utr_multiplier   = $param->{utr_multiplier};
  my $utr_extends      = $param->{utr_extends};

  my $flanks = {};
  my $utr_counts = {};
  $utr_counts->{3} = [];
  $utr_counts->{5} = [];

  open my $fh, '>', $filename;
  
  foreach my $transcript (@$transcripts) {
  
    (
      my $new_start, 
      my $new_end
    ) = compute_extended_transcript_coordinates({
      transcript       => $transcript, 
      utr_counts       => $utr_counts,
      unannotated_utrs => $unannotated_utrs,
      utr_multiplier   => $utr_multiplier,
      utr_extends      => $utr_extends,
      flanks           => $flanks
    });
    
    my $line_of_bed = join("\t", 
      $transcript->seq_region_name, 
      $new_start, 
      $new_end, 
      $transcript->stable_id
    );
  
    $fh->print($line_of_bed);
    $fh->print("\n");
  }
  $fh->close;

  return $flanks, $utr_counts;
}

sub compute_extended_transcript_coordinates {

  my $param = shift;
  lock_hash(%$param);

  my $transcript       = $param->{transcript};
  my $utr_counts       = $param->{utr_counts};
  my $unannotated_utrs = $param->{unannotated_utrs};
  my $utr_multiplier   = $param->{utr_multiplier};
  my $utr_extends      = $param->{utr_extends};
  my $flanks           = $param->{flanks};
  
  my $transcript_sid = $transcript->stable_id;

  # Computing UTR extension length
  $flanks->{$transcript_sid} = {};
  foreach my $end (5, 3) {
  
    my $extension_length = compute_extension_length({
      end            => $end, 
      transcript     => $transcript, 
      utr_counts     => $utr_counts, 
      utr_multiplier => $utr_multiplier,
      utr_extends    => $utr_extends,
    });
    $flanks->{$transcript_sid}{$end} = $extension_length;
  }

  my $new_start;
  my $new_end;
  # Check whether the transcript is on the postive/undefined or negative strand
  if ($transcript->strand() >= 0) {
    $new_start = $transcript->seq_region_start - $flanks->{$transcript_sid}{5};
    $new_end   = $transcript->seq_region_end   + $flanks->{$transcript_sid}{3};
  } else {
    $new_start = $transcript->seq_region_start - $flanks->{$transcript_sid}{3};
    $new_end   = $transcript->seq_region_end   + $flanks->{$transcript_sid}{5};
  }
  if ($new_start < 0) {
    $new_start = 0;
  }
  return ($new_start, $new_end)
}

=head2 compute_extension_length

  Arg1: "5" or "3"
  Arg2: Bio::EnsEMBL::Transcript object
  Arg3: Hash ref for UTR stats collection
  Arg4: Hash ref to collect stats
  Returntype: Length of extension in basepairs

=cut
sub compute_extension_length {

  my $param = shift;
  lock_hash(%$param);

  my $end            = $param->{end};
  my $transcript     = $param->{transcript};
  my $utr_counts     = $param->{utr_counts};
  my $utr_multiplier = $param->{utr_multiplier};
  my $utr_extends    = $param->{utr_extends};
  
  if($utr_extends->{$end} || $utr_multiplier || $unannotated_utrs->{$end}) {
    my $method = ($end == 5) ? 'five' : 'three';
    $method .= '_prime_utr';
    my $utr = $transcript->$method;

    if(defined $utr) {
      push @{$utr_counts->{$end}}, $utr->length;
      if(defined $utr_extends->{$end}) {
        return $utr_extends->{$end};
      }
      else{
        return $utr->length * $utr_multiplier;
      }
    } elsif (defined $unannotated_utrs->{$end}) {
      return $unannotated_utrs->{$end};
    }
  } else {
    return 0;
  }
}

