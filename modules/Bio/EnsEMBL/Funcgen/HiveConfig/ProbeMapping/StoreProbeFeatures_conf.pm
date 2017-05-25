package Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::StoreProbeFeatures_conf;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;           # Allow this particular config to use conditional dataflow and INPUT_PLUS
use base ('Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::Base');

my $max_allowed_hits_per_probe = 100;

sub pipeline_analyses {
    my $self = shift;
    
    return [
      {
          -logic_name  => 'start_store_probe_feature_chunk',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -flow_into => {
              'MAIN->A' => 'parse_exonerate',
              'A->MAIN' => 'done_store_probe_feature_chunk',
          },
      },
      {   -logic_name  => 'parse_exonerate',
          -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
              cmd => '
                import_parse_exonerate.pl \
                  --exonerate_file #chunk_name#_#type#.exonerate.txt \
                  --max_allowed_mismatches_per_hit 0 \
                  > #chunk_name#_#type#.exonerate_parsed.txt
              '
          },
          -flow_into => {
              MAIN => 'filter_features_from_promiscuous_probes',
          },
      },
      {   -logic_name  => 'filter_features_from_promiscuous_probes',
          -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -batch_size  => 50,
          -parameters => {
              cmd => '
                import_create_probe_feature_objects.pl \
                  --parsed_data #chunk_name#_#type#.exonerate_parsed.txt \
                  --promiscuous_hits #chunk_name#_#type#.promiscuous_hits.txt \
                  --accepted_hits #chunk_name#_#type#.probe_features.txt \
                  --max_allowed_hits_per_probe ' . $max_allowed_hits_per_probe
          },
          -flow_into => {
              MAIN => WHEN(
                  '#type# eq "genomic"'    => { 'store_probe_feature_objects'       => INPUT_PLUS },
                  '#type# eq "genomic"'    => { 'store_unmapped_objects'            => INPUT_PLUS },
                  '#type# eq "transcript"' => { 'project_transcript_hits_to_genome' => INPUT_PLUS },
              ),
          },
      },
      {   -logic_name  => 'project_transcript_hits_to_genome',
          -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -analysis_capacity => 60,
          -parameters => {
              cmd => '
                project_transcript_hits_to_genome.pl \
                  --probe_features #chunk_name#_#type#.probe_features.txt \
                  --registry #reg_conf# \
                  --species #species# \
                  --output_file #chunk_name#_#type#.probe_features_projected_from_transcript_coords.txt
              '
          },
          -flow_into => {
              MAIN => 'store_probe_feature_objects',
          },
      },
      {   -logic_name        => 'store_unmapped_objects',
          -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -priority          => 10,
          -analysis_capacity => 70,
          -parameters => {
              cmd => '
                store_unmapped_objects.pl \
                  --promiscuous_hits    #chunk_name#_#type#.promiscuous_hits.txt \
                  --registry            #reg_conf# \
                  --species             #species# \
                  --analysis_logic_name ProbeAlign_#type# \
                  --target_type         #type#
              '
          },
      },
      {   -logic_name        => 'store_probe_feature_objects',
          -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -priority          => 10,
          -analysis_capacity => 70,
          -parameters => {
              cmd => '
                store_probe_feature_objects.pl \
                  --probe_features #probe_feature_file# \
                  --registry #reg_conf# \
                  --species #species# \
                  --analysis_logic_name ProbeAlign_#type# \
                  --target_type #type#
              ',
              probe_feature_file => '#expr(#type# eq "genomic" ? #chunk_name#."_".#type#.".probe_features.txt" : #chunk_name#."_".#type#.".probe_features_projected_from_transcript_coords.txt" )expr#',
          },
      },
      {
          -logic_name  => 'done_store_probe_feature_chunk',
          -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      },
    ];
}

1;
