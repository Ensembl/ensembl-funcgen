#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Getopt::Long;
use Hash::Util qw( lock_hash );

=head1

associate_probes_to_transcripts.pl \
  --registry /homes/mnuhn/work_dir_probemapping/lib/ensembl-funcgen/registry.pm \
  --species homo_sapiens \
  --transcript_probe_features_overlaps /nfs/nobackup/ensembl/mnuhn/probe2transcript/transcript_probe_features_overlaps.bed \
  --flanks ./flanks.pl \
  --transcript_utr_file ./transcript_utr.pl


=cut

my $registry;
my $species;
my $transcript_file;

# Constants
my $debug = 0;
my $max_mismatches = 1;

GetOptions (
   'registry=s'                            => \$registry,
   'species=s'                             => \$species,
   'transcript_file=s'                     => \$transcript_file,
);

use Bio::EnsEMBL::Utils::Logger;
my $logger = Bio::EnsEMBL::Utils::Logger->new();

Bio::EnsEMBL::Registry->load_all($registry);

my $transcript_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'transcript');
my $transcripts        = $transcript_adaptor->fetch_all();

open my $transcript_probe_features_overlaps_fh, '<', $transcript_probe_features_overlaps;

my $whatsthis = associate_probes_to_transcripts({
  transcripts => $transcripts,
});

$transcript_probe_features_overlaps_fh->close;

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
  my $transcripts = $param->{transcripts};
  
  # Number of processed transcripts
  my $index = 0;
  # Hash ref with resuting overlaps
  my $object_transcript_hits = {};

  $logger->info('% Complete:', 0, 'append_date');
  foreach my $transcript (@$transcripts) {
    $last_pc = print_progress($index++, scalar @{$transcripts}, $last_pc);
    compute_hits(
      transcript              => $transcript,
      transcript_feature_info => $transcript_feature_info,
      object_transcript_hits  => $object_transcript_hits,
      unmapped_objects        => $unmapped_objects,
      unmapped_counts         => $unmapped_counts,
    );
  }
  return $object_transcript_hits;
}

=head2 print_progress

  Description: Prints a little blurb on stdout each % of the way through
  Arg1: Number of transcripts processd
  Arg2: Total number of transcripts
  Arg3: last percentage threshold to be printed out

=cut

sub print_progress {
  my ($index, $total, $last_pc) = @_;
  my $pc = int ((100 * $index) / $total);

  if ($pc > $last_pc) {
    $logger->info("$pc ", 0, 'append_date');
  }
  return $pc;
}

=head2

  Description: Converts hash ref linking objects to hits into a more conveninent
  hash ref object -> transcript -> array ref -> 2 integers
  Arg1: Bio::EnsEMBL::Transcript object
  Arg2: hash ref linking objects to hits
  Arg3: hash ref to be filled
  Arg4: array ref of Bio::EnsEMBL::UnmappedObject objects
  Arg5: hash ref with stats on unmapped objects
  Arg6: hash ref of global constants
  Arg7: File handle for detailed output

=cut

sub compute_hits {
  
  my $param = shift;
  lock_hash(%$param);

  my $transcript              = $param->{transcript};
  my $transcript_feature_info = $param->{transcript_feature_info};
  my $object_transcript_hits  = $param->{object_transcript_hits};
  my $unmapped_objects        = $param->{unmapped_objects};
  my $unmapped_counts         = $param->{unmapped_counts};
  
  my $compute_hits_2_parameters = {
      transcript              => $transcript,
      transcript_feature_info => $transcript_feature_info,
      object_transcript_hits  => $object_transcript_hits,
      unmapped_objects        => $unmapped_objects,
      unmapped_counts         => $unmapped_counts,
   };

  foreach my $object_id (keys %$transcript_feature_info) {
  
    $compute_hits_2_parameters->{object_id} = $object_id;  
    compute_hits_2($compute_hits_2_parameters);
  }
  return;
}

=head2 compute_hits_2

  Description: Converts hash ref linking an object to hits into a more conveninent
  hash ref object -> transcript -> array ref -> 2 integers
  Arg1: Ensembl DB ID
  Arg2: Bio::EnsEMBL::Transcript object
  Arg3: hash ref linking objects to hits
  Arg4: hash ref to be filled
  Arg5: array ref of Bio::EnsEMBL::UnmappedObject objects
  Arg6: hash ref with stats on unmapped objects
  Arg7: hash ref of global constants
  Arg8: File handle for detailed output

=cut

sub compute_hits_2 {

  my $param = shift;
  lock_hash(%$param);

  my $object_id               = $param->{object_id};
  my $transcript              = $param->{transcript};
  my $transcript_feature_info = $param->{transcript_feature_info};
  my $object_transcript_hits  = $param->{object_transcript_hits};
  my $unmapped_objects        = $param->{unmapped_objects};
  my $unmapped_object_counts  = $param->{unmapped_object_counts};

  # Count hits between object and transcript
  my $hits;
  if($options->{array_config}{probeset_arrays}) {
    $hits = scalar(keys %{$transcript_feature_info->{$object_id}});
  } else{
    $hits = scalar(@{$transcript_feature_info->{$object_id}});
  }

  my $probeset_size = $options->{probeset_sizes}{$object_id};

  # Verifiy whether the number of hits is above threshold
  if (
    ($options->{xref_object} eq 'ProbeSet' && ($hits / $probeset_size) >= $options->{mapping_threshold})
    || 
    ($hits && ($options->{xref_object} eq 'Probe'))
  ) {
    # Success: store the hit information
    store_sufficient_hit(
      $object_id               => $object_id,
      $transcript              => $transcript,
      $hits                    => $hits,
      $transcript_feature_info => $transcript_feature_info,
      $object_transcript_hits  => $object_transcript_hits,
    );
  } else {
    # Failure, store unmapped object information
    store_insufficient_hit(
      $object_id              => $object_id,
      $transcript             => $transcript,
      $hits                   => $hits,
      $probeset_size          => $probeset_size,
      $unmapped_objects       => $unmapped_objects,
      $unmapped_object_counts => $unmapped_object_counts,
    );
  }

  return;
}

=head2 store_sufficeient_hit

  Description: Converts hash ref linking an object to a transcript into a more conveninent
  hash ref object -> transcript -> array ref -> 2 integers
  Arg1: Ensembl DB ID
  Arg2: Bio::EnsEMBL::Transcript object
  Arg3: Number of hits
  Arg4: hash ref linking objects to hits
  Arg5: hash ref to be filled
   made up of probesets
  Arg6: hash ref of global constants
  Arg7: File handle for detailed output

=cut

sub store_sufficient_hit {

  my $param = shift;
  lock_hash(%$param);

  my $object_id               = $param->{object_id};
  my $transcript              = $param->{transcript};
  my $hits                    = $param->{hits};
  my $transcript_feature_info = $param->{transcript_feature_info};
  my $object_transcript_hits  = $param->{object_transcript_hits};

  # Count mismatches
  my $num_mismatch_hits = 0;
  if($options->{array_config}{probeset_arrays}) {
    foreach my $value (values %{$transcript_feature_info->{$object_id}}) {
      if ($value->[1] == $value->[0]) {
        $num_mismatch_hits += 1;
      }
    }
  } else{
    foreach my $value (@{$transcript_feature_info->{$object_id}}) {
      if ($value->[1]) {
        $num_mismatch_hits += 1;
      }
    }
  }

  $object_transcript_hits->{$object_id}{$transcript->stable_id} = [$hits, $num_mismatch_hits];
  if ($debug) {
    $logger->info("Sufficient hit\t$object_id\t".$transcript->stable_id()."\t$hits\t$num_mismatch_hits\n");
  }
  if ($options->{xref_object} eq 'Probe' && $hits == 1) {
    push @{$object_transcript_hits->{$object_id}{$transcript->stable_id}}, $transcript_feature_info->{$object_id}[0][0];
  }
}

=head2 store_insufficient_hit

  Description: Records a missed connection between object and transcript
  Arg1: Ensembl DB ID
  Arg2: Bio::EnsEMBL::Transcript object
  Arg3: Number of hits
  Arg4: Size of probeset
  Arg5: array ref of Bio::EnsEMBL::UnmappedObject objects
  Arg6: hash ref with stats on unmapped objects
   made up of probesets
  Arg7: hash ref of global constants
  Arg8: File handle for detailed output

=cut

sub store_insufficient_hit {

  my $param = shift;
  lock_hash(%$param);
  
  my $object_id              = $param->{object_id};
  my $transcript             = $param->{transcript};
  my $hits                   = $param->{hits};
  my $probeset_size          = $param->{probeset_size};
  my $unmapped_objects       = $param->{unmapped_objects};
  my $unmapped_object_counts = $param->{unmapped_object_counts};

  my $id_names = $object_id.'('.join(',', @{$options->{object_names}{$object_id}}).')';

  if ($debug) {
    $logger->info("$id_names\t".$transcript->stable_id."\tinsufficient\t$hits/$probeset_size in ProbeSet\n");
  }

  cache_and_load_unmapped_objects({
      identifier  => $transcript->stable_id,
      object_type => 'ProbeSet',
      object_id   => $object_id,
      summary     => 'Insufficient hits',
      description => "Insufficient number of hits $hits/$probeset_size in ProbeSet"
  });
}

=head2 cache_and_load_unmapped_objects

  Description: Stores an objects into the unmapped_object table, with some caching
  Arg1: hash ref of global constants
  Arg2: hash ref with stats on unmapped objects
  Arg3: array ref of Bio::EnsEMBL::UnmappedObject objects
  Arg4: stable id to potential hit
  Arg5: "Probe", "ProbeSet" or "ProbeFeature"
  Arg6: Ensembl ID of the object
  Arg7: short string description
  Arg8: longer string description

=cut

sub cache_and_load_unmapped_objects {

  my $param = shift;
  lock_hash(%$param);

  my $identifier       = $param->{identifier};
  my $object_type      = $param->{object_type};
  my $object_id        = $param->{object_id};
  my $summary          = $param->{summary};
  my $description      = $param->{description};

  my $um_obj = Bio::EnsEMBL::UnmappedObject->new(
    -type                => 'probe2transcript',
    -analysis            => $options->{analysis},
    -identifier          => $identifier,
    -summary             => $summary,
    -full_desc           => $description,
    -ensembl_object_type => $object_type,
    -ensembl_id          => $object_id,
#     -external_db_id      => $options->{transc_edb_id},
  );
  
  # Todo: Write this to file
  print "Todo: Write this to file\n";
  return;
}

