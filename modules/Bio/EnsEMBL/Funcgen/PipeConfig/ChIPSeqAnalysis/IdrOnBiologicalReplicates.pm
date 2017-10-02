package Bio::EnsEMBL::Funcgen::PipeConfig::ChIPSeqAnalysis::IdrOnBiologicalReplicates;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Funcgen::PipeConfig::ChIPSeqAnalysis::Base';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub pipeline_analyses {
    my $self = shift;

    return [
        {   -logic_name  => 'idr_on_biological_replicates',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               MAIN => 'br_idr_compute_peak_calling_inputs',
            },
        },
        {   -logic_name  => 'br_idr_compute_peak_calling_inputs',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::BrIdrComputePeakCallingInputs',
            -flow_into   => {
               2 => 'br_idr_start_idr',
               '3->A' => 'br_idr_start_align_all_read_files',
               'A->4' => 'br_idr_run_idr_peak_calling',
            },
        },
        @{
            $self->generate_parallel_alignment_analyses({
                start  => 'br_idr_start_align_all_read_files',
                prefix => 'br_idr_',
                suffix => '_all_read_files',
            })
        },
        @{
            $self->generate_parallel_alignment_analyses({
                start     => 'br_idr_start_idr',
                prefix    => 'br_idr_',
                suffix    => '_biological_replicates',
                flow_into => 'br_idr_run_permissive_swembl',
            })
        },
        {   -logic_name  => 'br_idr_run_permissive_swembl',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::RunPermissiveSWEmbl',
        },

        {   -logic_name  => 'br_idr_run_idr_peak_calling',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
    ]
}

1;

