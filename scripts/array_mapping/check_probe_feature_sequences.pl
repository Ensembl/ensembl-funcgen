#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Getopt::Long;

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Maxdepth = 3;

=head1

  Script to check probe features. For probe features that are perfect matches 
  it checks that the genome sequence they are on is identical to the probe 
  sequence.

  check_probe_feature_sequences.pl \
    --registry /homes/mnuhn/work_dir_probemapping/lib/ensembl-funcgen/registry.pm \
    --species saccharomyces_cerevisiae \
    --logic_name ProbeAlign_transcript \
    --check_probe_features_with_nontrivial_cigar_lines 1 \
    --max_check 10000

=cut

my $registry;
my $species;
my $max_check;
my $logic_name;
my $check_probe_features_with_nontrivial_cigar_lines;

GetOptions (
   'registry=s'   => \$registry,
   'species=s'    => \$species,
   'max_check=s'  => \$max_check,
   'logic_name=s' => \$logic_name,
   'check_probe_features_with_nontrivial_cigar_lines=s' => \$check_probe_features_with_nontrivial_cigar_lines,
);

Bio::EnsEMBL::Registry->load_all($registry);

my $funcgen_adaptor       = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'Funcgen' );
my $probe_feature_adaptor = Bio::EnsEMBL::Registry->get_adaptor( $species, 'Funcgen', 'ProbeFeature' );

use Bio::EnsEMBL::Utils::SqlHelper;
my $helper = Bio::EnsEMBL::Utils::SqlHelper->new(
  -DB_CONNECTION => $funcgen_adaptor->dbc
);

my $total_checked = 0;
my $total_failed  = 0;
my $total_probe_features_considered = 0;

my $limit_clause = '';
if ($max_check) {
  $limit_clause = "limit $max_check";
}

my @condition;

if ($logic_name) {
  push @condition, "logic_name = \"$logic_name\"";
}

if ($check_probe_features_with_nontrivial_cigar_lines) {
  push @condition, "cigar_line like '%D%'";
} else {
  push @condition, "cigar_line not like '%D%'";
}

my $where_clause = 'where ' . join ' and ', @condition;

my $sql = "select probe_feature_id from probe_feature join analysis using (analysis_id) $where_clause $limit_clause";

print "$sql\n";

$helper->execute_no_return(
  -SQL      => $sql,
  -CALLBACK => sub {
    my $row = shift;
    my $probe_feature_id = $row->[0];
    print check_probe_feature_with_id($probe_feature_id);
    return;
  },
);

my $msg = <<MESSAGE

total_probe_features_considered = $total_probe_features_considered

total_checked = $total_checked
total_failed  = $total_failed

MESSAGE
;

if ($total_failed>0) {
  die($msg);
}
print $msg;

sub check_probe_feature_with_id {

  my $probe_feature_id = shift;
  my $probe_feature = $probe_feature_adaptor->fetch_by_dbID($probe_feature_id);
  return check_probe_feature_and_update_stats($probe_feature);
}

sub check_probe_feature_and_update_stats {

  my $probe_feature = shift;
  
  $total_probe_features_considered++;
  (
    my $has_passed,
    my $message,
  ) = check_probe_feature($probe_feature);
  $total_checked++;

  if (! $has_passed) {
    $total_failed++;
    return $message;
  }
  return undef;
}

sub check_probe_feature {

  my $probe_feature = shift;
  
  my $parsed_cigar_string;
  eval {
    $parsed_cigar_string = parse_cigar_string($probe_feature->cigar_string);
  };
  if ($@) {
    die "Can't parse cigar string ". $probe_feature->cigar_string . " of probe feature with id: " . $probe_feature->dbID . "!\n";
  }
  
#   print Dumper($parsed_cigar_string);
  
  my $probe_feature_slice  = $probe_feature->slice;
  my $probe_feature_start  = $probe_feature->start;
  my $probe_feature_strand = $probe_feature->strand;
  
  my @matched_sequence_parts;
  
  ALIGNMENT_PART: foreach my $current_alignment_part (@$parsed_cigar_string) {
  
    my $is_a_match_region = $current_alignment_part->{type} eq '=';
    
    if (! $is_a_match_region) {
      next ALIGNMENT_PART;
    }
    
    my $matched_sequence_part = $probe_feature_slice->subseq(
      $probe_feature_start + $current_alignment_part->{start_on_alignment},
      $probe_feature_start + $current_alignment_part->{end_on_alignment},
      $probe_feature_strand,
    );
    push @matched_sequence_parts, $matched_sequence_part;
  }
  
  if ($probe_feature_strand == -1) {
    @matched_sequence_parts = reverse @matched_sequence_parts;
  }
  
  my $matched_sequence = join '', @matched_sequence_parts;
  
  my $probe_sequence = $probe_feature->probe->get_ProbeSequence->sequence;
  
  my $has_passed = uc($matched_sequence) eq uc($probe_sequence);
  
  my $summarise_test = 0;
  
  my $must_build_summary = $summarise_test || !$has_passed;

  my $summary;
  if ($must_build_summary) {
    $summary = 
        "Has passed: $has_passed\n"
      . "Cigar string:" . $probe_feature->cigar_string . "\n"
      . "Sequence in genome: " . (join '..', @matched_sequence_parts) . "\n"
      . "Probe sequence:     $probe_sequence\n";
  }
  
  if ($summarise_test) {
    print "--------------------------------------------------------------------\n";
    print "$summary";
  }
  return $has_passed, $summary;
}

sub parse_cigar_string {
  my $cigar_string = shift;
  
  my $parsed_cigar_lines = parse_cigar_string_generic($cigar_string);
  
  my $start_on_alignment = 0;
  my $end_on_alignment   = 0;

  foreach my $current_parsed_cigar_line (@$parsed_cigar_lines) {
  
    my $length = $current_parsed_cigar_line->{length};

    $end_on_alignment += $length;

    $current_parsed_cigar_line->{start_on_alignment} = $start_on_alignment;
    $current_parsed_cigar_line->{end_on_alignment}   = -1 + $end_on_alignment;
    
    use Hash::Util qw( lock_hash );
    lock_hash(%$current_parsed_cigar_line);
    
    $start_on_alignment += $length;
  }
  return $parsed_cigar_lines;
}

sub parse_cigar_string_generic {

  my $cigar_string = shift;
  
  my @x = $cigar_string =~ /(\d+)(\D)/g;
  # $VAR1 = [
  #   '14',
  #   '=',
  #   '669',
  #   'D',
  #   '35',
  #   '=',
  #   '181',
  #   'D',
  #   '1',
  #   '='
  # ];
  
  use List::MoreUtils qw( natatime );
  my @parsed_cigar_line;
  
  my $iterator = natatime 2, @x;
  while (my @pair = $iterator->()) {
  
    my $length = $pair[0];
    my $type   = $pair[1];
    
    push @parsed_cigar_line, {
      length => $length,
      type   => $type,
    };
  }
  return \@parsed_cigar_line;
}
