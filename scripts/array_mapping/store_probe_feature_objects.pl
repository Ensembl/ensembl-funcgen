#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Getopt::Long;

=head1 store_probe_feature_objects.pl

perl scripts/array_mapping/store_probe_feature_objects.pl \
  --probe_features /nfs/nobackup/ensembl/mnuhn/array_mapping/temp/homo_sapiens/probe_chunks/probe_chunk_2.fasta_genomic.probe_features.txt \
  --registry /nfs/gns/homes/mnuhn/work_dir_probemapping/lib/ensembl-funcgen/registry.pm \
  --species homo_sapiens \
  --analysis_logic_name ProbeAlign_genomic
  
=cut

my $probe_features;
my $registry;
my $species;
my $analysis_logic_name;
my $target_type;

GetOptions (
   'probe_features=s'      => \$probe_features,
   'registry=s'            => \$registry,
   'species=s'             => \$species,
   'analysis_logic_name=s' => \$analysis_logic_name,
   'target_type=s'         => \$target_type,
);

if (! -e $probe_features) {
  die("Can't find probe file ${probe_features}!");
}

use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->load_all($registry);

my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');

my $probe_feature_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'probefeature');
my $analysis_adaptor      = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'analysis');
my $probe_adaptor         = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'probe');
my $slice_adaptor         = Bio::EnsEMBL::Registry->get_adaptor($species, 'core',    'slice');
my $transcript_adaptor    = Bio::EnsEMBL::Registry->get_adaptor($species, 'core',    'transcript');

use Bio::EnsEMBL::Analysis;
my $analysis = Bio::EnsEMBL::Analysis->new(-logic_name => $analysis_logic_name);

$Data::Dumper::Sortkeys = 1;
# $Data::Dumper::Maxdepth = 3;

my $process_data = sub {

  my $raw_probe_feature = shift;
  
  use Hash::Util qw( lock_keys );
  lock_keys( %$raw_probe_feature );
  
  my $slice;
  my $hit_id;
  my $probe_feature_source;
  
  if ($target_type eq 'transcript') {
    # The next statement has to go
    return if ($raw_probe_feature eq " Only one genomic block!");
    my $transcript_stable_id = $raw_probe_feature->{t_id};
    my $transcript = $transcript_adaptor->fetch_by_stable_id($transcript_stable_id);
    
    if (! defined $transcript) {
      die("Can't find transcript for $transcript_stable_id");
    }
    $slice  = $transcript->slice;
    $hit_id = $transcript->stable_id;
    $probe_feature_source = 'transcript';
    
  } else {
    $slice = $slice_adaptor->fetch_by_name($raw_probe_feature->{t_id});
    $hit_id = $slice->seq_region_name;
    $probe_feature_source = 'genomic';
  }
  
  my $probe_with_this_sequence = $probe_adaptor->fetch_all_by_probe_sequence_id($raw_probe_feature->{probe_seq_id});
  
  foreach my $current_probe (@$probe_with_this_sequence) {
  
    use Bio::EnsEMBL::Funcgen::ProbeFeature;
    my $probe_feature = Bio::EnsEMBL::Funcgen::ProbeFeature->new(
      -PROBE         => $current_probe,
      -MISMATCHCOUNT => $raw_probe_feature->{total_mismatches},
      -START         => $raw_probe_feature->{t_start},
      -END           => $raw_probe_feature->{t_end},
      -STRAND        => $raw_probe_feature->{t_strand},
      -ANALYSIS      => $analysis,
      -CIGAR_STRING  => $raw_probe_feature->{cigar_line},
      -slice         => $slice,
      -hit_id        => $hit_id,
      -source        => $probe_feature_source,
    );

    $probe_feature_adaptor->store($probe_feature);
  }
};

use Bio::EnsEMBL::Funcgen::Parsers::DataDumper;
my $parser = Bio::EnsEMBL::Funcgen::Parsers::DataDumper->new;

# This is what the parser will return. If it is not "used", the inheritance 
# in Bio::EnsEMBL::Funcgen::Probe to Bio::EnsEMBL::Storable will not work.
#
use Bio::EnsEMBL::Funcgen::Probe;

$parser->parse({
  data_dumper_file => $probe_features,
  call_back        => $process_data,
});
