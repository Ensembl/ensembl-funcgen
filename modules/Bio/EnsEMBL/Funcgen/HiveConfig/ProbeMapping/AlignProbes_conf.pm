package Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::AlignProbes_conf;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;           # Allow this particular config to use conditional dataflow and INPUT_PLUS
use base ('Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::Base');

sub pipeline_analyses {
    my $self = shift;
    
    return [
      {
          -logic_name  => 'start_align_probes',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -flow_into => {
              MAIN => 'job_factory_probe_align',
          },
      },
        {   -logic_name  => 'job_factory_probe_align',
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
                2 => [ 
                  { 'multiply_probe_align_jobs' => INPUT_PLUS },
                ],
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
                MAIN => 'start_store_probe_features',
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
                MAIN => 'start_store_probe_features',
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
                MAIN => 'start_store_probe_features',
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
                MAIN => 'start_store_probe_features',
            },
        },
      {
          -logic_name  => 'start_store_probe_features',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      },
    ];
}

1;
