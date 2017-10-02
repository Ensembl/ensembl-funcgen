package Bio::EnsEMBL::Funcgen::PipeConfig::ChIPSeqAnalysis::IdrOnTechnicalReplicates;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Funcgen::PipeConfig::ChIPSeqAnalysis::Base';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub pipeline_analyses {
    my $self = shift;

    return [
        {   -logic_name  => 'idr_on_technical_replicates',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               MAIN => 'tr_idr_compute_peak_calling_inputs',
            },
        },
        {   -logic_name  => 'tr_idr_compute_peak_calling_inputs',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::TrIdrComputePeakCallingInputs',
            -flow_into   => {
               2 => 'tr_idr_start_idr',
               '3->A' => 'tr_idr_start_align_all_read_files',
               'A->4' => 'tr_idr_run_idr',
            },
        },
        @{
            $self->generate_parallel_alignment_analyses({
                start  => 'tr_idr_start_align_all_read_files',
                prefix => 'tr_idr_',
                suffix => '_all_read_files',
            })
        },
        @{
            $self->generate_parallel_alignment_analyses({
                start  => 'tr_idr_start_idr',
                prefix => 'tr_idr_',
                suffix => '_technical_replicates',
                flow_into => 'tr_idr_run_permissive_swembl',
            })
        },
        {   -logic_name  => 'tr_idr_run_permissive_swembl',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::RunPermissiveSWEmbl',
        },
        {   -logic_name  => 'tr_idr_run_idr',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
    ]
}

1;

