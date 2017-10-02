package Bio::EnsEMBL::Funcgen::PipeConfig::ChIPSeqAnalysis::NoIdr;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Funcgen::PipeConfig::ChIPSeqAnalysis::Base';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub pipeline_analyses {
    my $self = shift;

    return [
        {   -logic_name  => 'no_idr',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               MAIN => 'no_idr_compute_peak_calling_inputs',
            },
        },
        {   -logic_name  => 'no_idr_compute_peak_calling_inputs',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::NoIdrComputePeakCallingInputs',
            -flow_into   => {
               '2->A' => 'no_idr_start_align_all_read_files',
               'A->4' => 'no_idr_call_peaks',
            },
        },
        @{
            $self->generate_parallel_alignment_analyses({
                start  => 'no_idr_start_align_all_read_files',
                prefix => 'no_idr_',
                suffix => '_all_read_files',
            })
        },
        {   -logic_name  => 'no_idr_call_peaks',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        {   -logic_name => 'no_idr_done',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
    ]
}

1;

