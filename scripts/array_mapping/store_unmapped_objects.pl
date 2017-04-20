#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Getopt::Long;

=head1 store_unmapped_objects.pl

perl scripts/array_mapping/store_unmapped_objects.pl \
  --promiscuous_hits /nfs/nobackup/ensembl/mnuhn/array_mapping/temp/homo_sapiens/probe_chunks/6/6/6/probe_chunk_1666.fasta_genomic.promiscuous_hits.txt \
  --registry /nfs/gns/homes/mnuhn/work_dir_probemapping/lib/ensembl-funcgen/registry.pm \
  --species homo_sapiens \
  --target_type genomic \
  --analysis_logic_name ProbeAlign_genomic
  
=cut

my $promiscuous_hits;
my $registry;
my $species;
my $analysis_logic_name;
my $target_type;
my $max_allowed_hits_per_probe = 100;

GetOptions (
   'promiscuous_hits=s'    => \$promiscuous_hits,
   'registry=s'            => \$registry,
   'species=s'             => \$species,
   'analysis_logic_name=s' => \$analysis_logic_name,
   'target_type=s'         => \$target_type,
   'max_allowed_hits_per_probe=s' => \$max_allowed_hits_per_probe,
);

if (! -e $promiscuous_hits) {
  die("Can't find file ${promiscuous_hits}!");
}

use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->load_all($registry);

my $analysis_adaptor        = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'analysis');
my $unmapped_object_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'unmappedobject');
my $probe_adaptor           = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'probe');

use Bio::EnsEMBL::Analysis;
my $analysis = Bio::EnsEMBL::Analysis->new(-logic_name => $analysis_logic_name);

my $process_data = sub {

  my $count_probe_match_counts = shift;
  
  foreach my $current_probe_seq_id (keys %$count_probe_match_counts) {
    
    my $probe_with_this_sequence = $probe_adaptor->fetch_all_by_probe_sequence_id($current_probe_seq_id);
    
    foreach my $current_probe (@$probe_with_this_sequence) {
      my $match_count = $count_probe_match_counts->{$current_probe_seq_id};

      if ($match_count>$max_allowed_hits_per_probe) {
      
        use Bio::EnsEMBL::UnmappedObject;
        my $unmapped_object = Bio::EnsEMBL::UnmappedObject->new (
          -type                => 'array_mapping',
          -analysis            => $analysis,
          -ensembl_id          => $current_probe->dbID,
          -ensembl_object_type => 'Probe',
          -external_db_id      => undef,
          -identifier          => $target_type,
          -summary             => 'Promiscuous probe',
          -full_desc           => "Probe exceeded maximum allowed number of mappings. ($match_count matches)"
        );
        $unmapped_object_adaptor->store($unmapped_object);
      }
    }
  }
};

use Bio::EnsEMBL::Funcgen::Parsers::DataDumper;
my $parser = Bio::EnsEMBL::Funcgen::Parsers::DataDumper->new;

$parser->parse({
  data_dumper_file => $promiscuous_hits,
  call_back        => $process_data,
});
