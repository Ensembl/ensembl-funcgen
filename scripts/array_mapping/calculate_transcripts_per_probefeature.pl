#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Getopt::Long;

=head1

date
calculate_transcripts_per_probefeature.pl \
  --registry /homes/mnuhn/work_dir_probemapping/lib/ensembl-funcgen/registry.pm \
  --species homo_sapiens \
  --transcripts_per_probefeature_file ./transcripts_per_probefeature.pl
date

=cut

my $debug;

my $registry;
my $species;
my $transcripts_per_probefeature_file;

GetOptions (
   'registry=s'                     => \$registry,
   'species=s'                      => \$species,
   'transcripts_per_probefeature_file=s' => \$transcripts_per_probefeature_file,
);

use Bio::EnsEMBL::Utils::Logger;
my $logger = Bio::EnsEMBL::Utils::Logger->new();

Bio::EnsEMBL::Registry->load_all($registry);

my $funcgen_db_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen', 'transcript');
my $transcripts_per_probefeature = get_transcripts_per_probefeature($funcgen_db_adaptor);

open my $out, '>', $transcripts_per_probefeature_file;
$out->print(Dumper($transcripts_per_probefeature));
$out->close;

=head2 get_transcripts_per_probefeature

  Description: counts how many transcripts overlap a probe feature from the xref table
  Arg1: Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor object
  Arg2: Hash ref containing:
  - array_names: array ref
  - transc_edb_id: id of the core database in the external_db table of the xref system
  Arg3: Output file handle

=cut

sub get_transcripts_per_probefeature {
  my $xref_db = shift;
  my $sql = 'select probe_feature_id, stable_id from probe_feature_transcript';
  
  my $sth = $xref_db->dbc->prepare($sql);
  $sth->{mysql_use_result}=1; 
  $sth->execute();

  my %transcripts_per_probefeature;
  if ($debug) {
    $logger->info($sql."\n");
  }
  
  PROBE_FEATURE_ID:
  while(my $hash = $sth->fetchrow_hashref) {
    
    my $probe_feature_id = $hash->{probe_feature_id};
    my $transcript_sid   = $hash->{stable_id};
    
    if ($debug) {
      $logger->info("ALIGNMENT\t$probe_feature_id\t$transcript_sid\n");
    }

    if (
         (! exists $transcripts_per_probefeature{$probe_feature_id})
      || (! exists $transcripts_per_probefeature{$probe_feature_id}{$transcript_sid})
    ) {
      $transcripts_per_probefeature{$probe_feature_id} = {
        $transcript_sid => 1
      };
      next PROBE_FEATURE_ID;
    }
    $transcripts_per_probefeature{$probe_feature_id}{$transcript_sid}++;
  }

  return \%transcripts_per_probefeature;
}
