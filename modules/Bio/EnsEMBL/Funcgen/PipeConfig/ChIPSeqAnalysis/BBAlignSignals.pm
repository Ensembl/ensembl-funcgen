package Bio::EnsEMBL::Funcgen::PipeConfig::ChIPSeqAnalysis::BBAlignSignals;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Funcgen::PipeConfig::ChIPSeqAnalysis::Base';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub pipeline_analyses {
    my $self = shift;
    return [
        {
            -logic_name  => 'start_align_signals',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               MAIN => 'seed_signal_alignments',
            },
        },
        {   -logic_name  => 'seed_signal_alignments',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::SeedSignalAlignments',
            -flow_into   => {
               2 => 'start_aligning_signals',
            },
        },
        @{
            $self->generate_parallel_alignment_analyses({
                start  => 'start_aligning_signals',
                prefix => '',
                suffix => '_signals',
            })
        },
    ];
}

1;

