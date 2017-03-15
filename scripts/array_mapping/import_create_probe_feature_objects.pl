#!/usr/bin/env perl

use strict;
use Data::Dumper;
use Getopt::Long;
use Hash::Util qw( lock_keys );

=head1

perl scripts/array_mapping/import_create_probe_feature_objects.pl \
  --parsed_data /nfs/nobackup/ensembl/mnuhn/array_mapping/temp/homo_sapiens/probe_chunks/probe_chunk_2.fasta_genomic.exonerate_parsed.txt \
  --promiscuous_hits /nfs/nobackup/ensembl/mnuhn/array_mapping/temp/homo_sapiens/probe_chunks/probe_chunk_2.fasta_genomic.promiscuous_hits.txt \
  --accepted_hits /nfs/nobackup/ensembl/mnuhn/array_mapping/temp/homo_sapiens/probe_chunks/probe_chunk_2.fasta_genomic.accepted_hits.txt \
  --max_allowed_hits_per_probe 100

=cut

my $parsed_data;
my $output_file;
my $max_allowed_hits_per_probe;
my $promiscuous_hits;
my $accepted_hits;

GetOptions (
   'parsed_data=s'                => \$parsed_data,
   'promiscuous_hits=s'           => \$promiscuous_hits,
   'accepted_hits=s'              => \$accepted_hits,
   'max_allowed_hits_per_probe=s' => \$max_allowed_hits_per_probe,
);

open my $promiscuous_hits_fh, '>', $promiscuous_hits;
open my $accepted_hits_fh,    '>', $accepted_hits;

my %count_probe_match_counts;

# $count_probe_match_counts{'foo'} = 'bar';

my $count_probe_matches = sub {

  my $probe_hit_list = shift;
  
  foreach my $current_probe_hit (@$probe_hit_list) {
  
    my $probe_seq_id = $current_probe_hit->{probe_seq_id};
  
    if (!exists $count_probe_match_counts{$probe_seq_id}) {
      $count_probe_match_counts{$probe_seq_id} = 1;
    } else {
      $count_probe_match_counts{$probe_seq_id}++
    }
  }
};

my $process_probe_data = sub {

  my $probe_hit_list = shift;
  
  foreach my $current_probe_hit (@$probe_hit_list) {
  
    lock_keys( %$current_probe_hit);
  
    my $probe_seq_id = $current_probe_hit->{probe_seq_id};
    
    if (!exists $count_probe_match_counts{$probe_seq_id}) {
      die;
    }
    
    my $number_of_probe_matches = $count_probe_match_counts{$probe_seq_id};
    
    if ($number_of_probe_matches<$max_allowed_hits_per_probe) {
      $accepted_hits_fh->print(Dumper($current_probe_hit));
    }
  }
};

use Bio::EnsEMBL::Funcgen::Parsers::DataDumper;
my $parser = Bio::EnsEMBL::Funcgen::Parsers::DataDumper->new;

$parser->parse({
  data_dumper_file => $parsed_data,
  call_back        => $count_probe_matches,
});

$promiscuous_hits_fh->print(Dumper(\%count_probe_match_counts));

$parser->parse({
  data_dumper_file => $parsed_data,
  call_back        => $process_probe_data,
});
