package Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::AlignProbes_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::Base');

sub pipeline_analyses {
    my $self = shift;
    
    return [
      {
          -logic_name  => 'start',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -input_ids => [
            { species => 'homo_sapiens', },
          ],
          -flow_into => {
              MAIN => 'job_factory_probe_align',
          },
      },
        {   -logic_name  => 'job_factory_probe_align',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::FastaFactory',
            -parameters => {
              inputfile        => '#unmapped_sequences_file#',
              output_dir       => '#tempdir#/#species#/probe_chunks',
              max_chunk_length => 1000,
              output_prefix    => 'probe_chunk_',
              output_suffix    => '.fasta',
              hash_directories => 1,
            },
            -flow_into => {
                MAIN => [ 
                  'probe_align_genomic',
                  'probe_align_transcript', 
                ],
            },
        },
        {   -logic_name  => 'probe_align_genomic',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::ProbeAlign',
            -rc_name     => '8Gb_job',
            -analysis_capacity => 20,
            -batch_size        => 100,
            -parameters => {
                QUERYSEQS    => '#chunk_name#',
                TARGETSEQS   => '#toplevel_sequences_file#',
                mapping_type => 'genomic', # genomic or transcript
            },
            -flow_into => {
                MAIN => 'probe_align_genomic_himem',
            },
        },
        {   -logic_name  => 'probe_align_genomic_himem',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::ProbeAlign',
            -rc_name     => '64Gb_job',
            -analysis_capacity => 20,
            -batch_size        => 100,
            -parameters => {
                QUERYSEQS    => '#chunk_name#',
                TARGETSEQS   => '#toplevel_sequences_file#',
                mapping_type => 'genomic', # genomic or transcript
            },
        },
        {   -logic_name  => 'probe_align_transcript',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::ProbeAlign',
            -batch_size => 100,
            -rc_name     => '8Gb_job',
            -parameters => {
                QUERYSEQS    => '#chunk_name#',
                TARGETSEQS   => '#gene_sequences_file#',
                mapping_type => 'transcript', # genomic or transcript
            },
            -flow_into => {
                MAIN => 'probe_align_transcript_himem',
            },
        },
        {   -logic_name  => 'probe_align_transcript_himem',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::ProbeAlign',
            -rc_name     => '64Gb_job',
            -batch_size => 100,
            -parameters => {
                QUERYSEQS    => '#chunk_name#',
                TARGETSEQS   => '#gene_sequences_file#',
                mapping_type => 'transcript', # genomic or transcript
            },
        },
    ];
}

1;
