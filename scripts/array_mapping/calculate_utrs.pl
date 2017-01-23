#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Getopt::Long;

=head1

date
calculate_utrs.pl \
  --registry /homes/mnuhn/work_dir_probemapping/lib/ensembl-funcgen/registry.pm \
  --species homo_sapiens \
  --transcript_utr_file ./transcript_utr2.pl \
date

=cut

my $registry;
my $species;
my $transcript_utr_file;

GetOptions (
   'registry=s'            => \$registry,
   'species=s'             => \$species,
   'transcript_utr_file=s' => \$transcript_utr_file,
);

use Bio::EnsEMBL::Utils::Logger;
my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

Bio::EnsEMBL::Registry->load_all($registry);

my $transcript_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'transcript');

$logger->info("Fetching all transcripts for $species\n");
my $transcripts        = $transcript_adaptor->fetch_all();
$logger->info("Done.\n");
$logger->info("Got ". scalar @$transcripts ." transcripts.\n");

$logger->info("Calculating utrs.\n");
my $unannotated_utrs = calculate_utrs($transcripts);
$logger->info("Done.\n");

$logger->info("Storing results in $transcript_utr_file.\n");
open my $out, '>', $transcript_utr_file;
$out->print(Dumper($unannotated_utrs));
$out->close;
$logger->info("Done.\n");

$logger->finish_log;

=head2 calculate_utrs

  Description: Computes the extension required on each gene. Fills in unannotated_utrs hashref
  Arg1: hashref, possibly with integer values associated to "5" and "3" constants.
  Arg2: Array ref of Bio::EnsEMBL::Transcript objects

=cut
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw (median mean);
sub calculate_utrs {
  my $transcripts = shift;
  
  my $unannotated_utrs;
  my %lengths = (
    5 => [],
    3 => []
  );

  # Get a list of all utr lengths from all transcripts, where defined
  foreach my $transcript(@$transcripts) {
    if(! defined $unannotated_utrs->{3} && $transcript->three_prime_utr) {
      push @{$lengths{3}}, $transcript->three_prime_utr->length;
    }
    if(! defined $unannotated_utrs->{5} && $transcript->five_prime_utr) {
      push @{$lengths{5}}, $transcript->five_prime_utr->length;
    }
  }

  # Compute a typical utr length from the measured values
  foreach my $side (5,3) {
    if(! defined $unannotated_utrs->{$side}) {
      my $count = scalar @{$lengths{$side}};
      my $zero_count = (scalar @$transcripts) - $count;
      $logger->info("Seen $count ${side}' UTRs, $zero_count have length 0\n");

      if($count) {
        my $mean   = round(mean($lengths{$side}));
        my $median = round(median($lengths{$side}));
        $unannotated_utrs->{$side}  = ($mean > $median)  ? $mean : $median;
        $logger->info("Calculated default unannotated ${side}' UTR length:\t$unannotated_utrs->{$side}\n");
      } else{
        use Carp;
        warn("Found no ${side}' UTRs, you must specify a -unannotated_${side}_utr");
        $unannotated_utrs->{$side}  = 0;
      }
    }
  }
  return $unannotated_utrs;
}

sub round {
  my $number = shift;

  my ($rounded, $remainder) = split/\./, $number;
  if ($remainder =~ /^[5-9]/) {
    $rounded++;
  }
  return $rounded;
}



