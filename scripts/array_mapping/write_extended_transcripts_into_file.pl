#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Getopt::Long;
use Hash::Util qw( lock_hash );

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
   'registry=s'                       => \$registry,
   'species=s'                        => \$species,
   'unannotated_utrs=s'               => \$unannotated_utrs_file,
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

