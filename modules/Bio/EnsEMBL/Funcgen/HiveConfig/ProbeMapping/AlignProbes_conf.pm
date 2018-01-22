package Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::AlignProbes_conf;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;           # Allow this particular config to use conditional dataflow and INPUT_PLUS
use base ('Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::Base');

my $max_allowed_hits_per_probe = 100;

sub pipeline_analyses {
    my $self = shift;
    
    return [
      {
          -logic_name  => 'start_align_probes',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -flow_into => {
              MAIN => 'truncate_probe_feature_tables',
          },
      },
        {
            -logic_name  => 'truncate_probe_feature_tables',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
            -parameters => {
                sql     => [
                  "truncate probe_feature;",
                  "truncate probe_transcript;",
                  "truncate probe_set_transcript;",
                  "truncate probe_feature_transcript;",
                  "truncate unmapped_object;",
                  "truncate unmapped_reason;",
                ],
                db_conn => 'funcgen:#species#',
            },
            -flow_into => {
               MAIN => 'switch_probe_feature_table_to_innodb',
            },
        },
        {
            -logic_name  => 'switch_probe_feature_table_to_innodb',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
            -parameters => {
                sql     => [
                  "ALTER TABLE probe_feature ENGINE=InnoDB;",
                ],
                db_conn => 'funcgen:#species#',
            },
            -flow_into => {
               'MAIN->A' => 'split_into_probe_chunks',
               'A->MAIN' => 'done_processing_probe_chunks',
            },
        },
        {   -logic_name  => 'split_into_probe_chunks',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::FastaFactory',
            -parameters => {
              inputfile        => '#unmapped_sequences_file#',
              output_dir       => '#tempdir#/#species#/probe_chunks',
              max_chunk_length => 100000,
              output_prefix    => 'probe_chunk_',
              output_suffix    => '.fasta',
              hash_directories => 1,
            },
            -flow_into => {
                2 => [ { 'multiply_probe_align_jobs' => INPUT_PLUS }, ],
#                 '2->A' => [ { 'multiply_probe_align_jobs' => INPUT_PLUS }, ],
#                 'A->2' => [ { 'done_processing_probe_chunks' => undef      }, ]
            },
        },
      {
          -logic_name  => 'multiply_probe_align_jobs',
          -module     => 'Bio::EnsEMBL::Funcgen::Hive::AddJobParameter',
          -parameters => {
            add => {
              species => '#species#',
            }
          },
          -flow_into => {
              2 => [
                'start_align_probes_genomic',
                'start_align_probes_transcript'
               ],
          },
      },
      {
          -logic_name  => 'start_align_probes_genomic',
          -module     => 'Bio::EnsEMBL::Funcgen::Hive::AddJobParameter',
          -parameters => {
            add => {
              type    => 'genomic',
            }
          },
          -flow_into => {
              2 => 'probe_align_genomic',
          },
      },
      {
          -logic_name  => 'start_align_probes_transcript',
          -module     => 'Bio::EnsEMBL::Funcgen::Hive::AddJobParameter',
          -parameters => {
            add => {
              type    => 'transcript',
            }
          },
          -flow_into => {
              2 => 'probe_align_transcript',
          },
      },
        {   -logic_name  => 'probe_align_genomic',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -rc_name     => '8Gb_job',
#             -analysis_capacity => 1,
            -batch_size        => 10,
            -parameters => {
                cmd => '
                  exonerate \
                  --querytype dna \
                  --targettype dna \
                  --query #chunk_name# \
                  --target #toplevel_sequences_file# \
                  --showsugar false --showvulgar false --showalignment false --ryo "RESULT: %S %pi %ql %tl %em scores:{%Ps:}\\n" \
                  > #chunk_name#_genomic.exonerate.txt
                '
            },
            -flow_into => {
                MAIN => 'start_store_probe_feature_chunk',
                '-1' => 'probe_align_genomic_himem',
            },
        },
        {   -logic_name  => 'probe_align_genomic_himem',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -rc_name     => '64Gb_job',
            -batch_size        => 10,
            -parameters => {
                cmd => '
                  exonerate \
                  --querytype dna \
                  --targettype dna \
                  --query #chunk_name# \
                  --target #toplevel_sequences_file# \
                  --showsugar false --showvulgar false --showalignment false --ryo "RESULT: %S %pi %ql %tl %em scores:{%Ps:}\\n" \
                  > #chunk_name#_genomic.exonerate.txt
                '
            },
            -flow_into => {
                MAIN => 'start_store_probe_feature_chunk',
            },
        },
        {   -logic_name  => 'probe_align_transcript',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -batch_size => 10,
            -rc_name     => '8Gb_job',
#             -analysis_capacity => 1,
            -parameters => {
                cmd => '
                  exonerate \
                  --querytype dna \
                  --targettype dna \
                  --query #chunk_name# \
                  --target #gene_sequences_file# \
                  --showsugar false --showvulgar false --showalignment false --ryo "RESULT: %S %pi %ql %tl %em scores:{%Ps:}\\n" \
                  > #chunk_name#_transcript.exonerate.txt
                '
            },
            -flow_into => {
                MAIN => 'start_store_probe_feature_chunk',
                '-1' => 'probe_align_transcript_himem',
            },
        },
        {   -logic_name  => 'probe_align_transcript_himem',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -batch_size => 10,
            -rc_name     => '64Gb_job',
            -parameters => {
                cmd => '
                  exonerate \
                  --querytype dna \
                  --targettype dna \
                  --query #chunk_name# \
                  --target #gene_sequences_file# \
                  --showsugar false --showvulgar false --showalignment false --ryo "RESULT: %S %pi %ql %tl %em scores:{%Ps:}\\n" \
                  > #chunk_name#_transcript.exonerate.txt
                '
            },
            -flow_into => {
                MAIN => 'start_store_probe_feature_chunk',
            },
        },
      {
          -logic_name  => 'start_store_probe_feature_chunk',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      },
      {
          -logic_name  => 'done_processing_probe_chunks',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -flow_into => {
              MAIN => 'delete_redundant_probe_features',
          },
      },
      {   -logic_name        => 'delete_redundant_probe_features',
          -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
              cmd => '
                delete_redundant_probe_features.pl \
                  --registry #reg_conf# \
                  --species #species#
              ',
          },
          -flow_into => {
              MAIN => 'delete_probe_features_on_transcripts_from_methylation_arrays',
          },
      },
      {
          -logic_name  => 'delete_probe_features_on_transcripts_from_methylation_arrays',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
          -parameters => {
              sql     => [
                'delete from probe_feature using probe_feature join probe using (probe_id) join array_chip using (array_chip_id) join array using (array_id) where array.format="METHYLATION" and probe_feature.source="transcript"',
              ],
              db_conn => 'funcgen:#species#',
          },
          -flow_into => {
              MAIN => 'delete_probe_features_from_known_promiscuous_probes',
          },
      },
      {
          -logic_name  => 'delete_probe_features_from_known_promiscuous_probes',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
          -parameters => {
              sql     => [
                'drop table if exists temp_promiscuous_probes_with_features;',

                "
                create table temp_promiscuous_probes_with_features as 
                  select
                    distinct unmapped_object.ensembl_id as probe_id
                  from
                    unmapped_object
                    join unmapped_reason on (
                      unmapped_object.unmapped_reason_id=unmapped_reason.unmapped_reason_id and summary_description = 'Promiscuous probe'
                    )
                    join probe_feature on (probe_id = ensembl_id and ensembl_object_type = 'Probe');
                ",
                
                'delete from probe_feature using probe_feature join temp_promiscuous_probes_with_features using(probe_id);',

                'drop table if exists temp_promiscuous_probes_with_features;',
              ],
              db_conn => 'funcgen:#species#',
          },
          -flow_into => {
              MAIN => 'delete_probe_features_from_promiscuous_probes',
          },
      },
      # This will typically pick up on probes that have scored less than the 
      # $max_allowed_hits_per_probe probe features for genomic hits and less
      # than $max_allowed_hits_per_probe for transcript hits. But when both
      # numbers are combined, they exceed $max_allowed_hits_per_probe and
      # the probe is classified as promiscuous after all.
      # 
      # This analysis will find probes of this kind, remove their probe 
      # features and create an unmapped object summarising the number of
      # probe features from transcript and genomic matches.
      #
      {   -logic_name        => 'delete_probe_features_from_promiscuous_probes',
          -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
              cmd => '
                delete_promiscuous_probe_features.pl \
                  --species  #species# \
                  --registry #reg_conf# \
                  --analysis_logic_name ProbeAlign_genomic \
                  --max_allowed_hits_per_probe ' . $max_allowed_hits_per_probe,
          },
      },
    ];
}

1;
