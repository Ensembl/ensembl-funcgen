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
    --max_check 10000

=cut

my $registry;
my $species;
my $max_check;
my $logic_name;

GetOptions (
   'registry=s'   => \$registry,
   'species=s'    => \$species,
   'max_check=s'  => \$max_check,
   'logic_name=s' => \$logic_name,
);

Bio::EnsEMBL::Registry->load_all($registry);

my $funcgen_adaptor       = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'Funcgen' );
my $probe_feature_adaptor = Bio::EnsEMBL::Registry->get_adaptor( $species, 'Funcgen', 'ProbeFeature' );

use Bio::EnsEMBL::Utils::SqlHelper;
my $helper = Bio::EnsEMBL::Utils::SqlHelper->new(
  -DB_CONNECTION => $funcgen_adaptor->dbc
);

my $total_checked = 0;
my $total_skipped = 0;
my $total_failed  = 0;
my $total_probe_features_considered = 0;

my $limit_clause = '';
if ($max_check) {
  $limit_clause = "limit $max_check";
}

my $where_clause = '';
if ($logic_name) {
  $where_clause = "where logic_name = \"$logic_name\"";
}

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
total_skipped = $total_skipped
total_failed  = $total_failed

MESSAGE
;

if ($total_failed>0) {
  die($msg);
}
print $msg;

sub check_probe_feature_with_id {

  my $probe_feature_id = shift;
  my $probe_feature_with_better_slice = $probe_feature_adaptor->fetch_by_dbID($probe_feature_id);
  return check_probe_feature($probe_feature_with_better_slice);
}

sub check_probe_feature {

  my $probe_feature = shift;
  
  $total_probe_features_considered++;

  # Either 'ProbeAlign_transcript' or 'ProbeAlign_genomic'
  #
  if ($probe_feature->analysis->logic_name eq 'ProbeAlign_transcript') {
  
    my $probe_feature_alignment_has_deletions = $probe_feature->cigar_string =~ /D/;
  
    # If there are deletions, then the sequences are not expected to be identical.
    #
    if ($probe_feature_alignment_has_deletions) {
      $total_skipped++;
      return;
    }
  }
  
  # The pipeline doesn't currently allow mismatches.
  #
  if ($probe_feature->mismatchcount != 0) {
    print Dumper($probe_feature);
    die;
  }
  my $matched_sequence = $probe_feature->slice->subseq(
    $probe_feature->start,
    $probe_feature->end,
    $probe_feature->strand,
  );

  my $probe_sequence = uc($probe_feature->probe->fetch_ProbeSequence->sequence);
  
  my $match_ok = $matched_sequence eq $probe_sequence;
  $total_checked++;
  
  if (! $match_ok) {
    $total_failed++;
    return "Not ok: " . $probe_feature->analysis->logic_name . " " . $probe_feature->dbID . " $matched_sequence !=  $probe_sequence " . $probe_feature->cigar_string . " " . $probe_feature->strand  . "\n";
  }
}
