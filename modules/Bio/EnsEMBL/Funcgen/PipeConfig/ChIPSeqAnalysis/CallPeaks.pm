package Bio::EnsEMBL::Funcgen::PipeConfig::ChIPSeqAnalysis::CallPeaks;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Funcgen::PipeConfig::ChIPSeqAnalysis::Base';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub pipeline_analyses {
    my $self = shift;
    return [
       {   -logic_name  => 'start_call_peaks',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               'MAIN->A' => 'call_peaks',
               'A->MAIN' => 'done_call_peaks',
            },
        },
        {   -logic_name  => 'call_peaks',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::CallPeaks',
            -flow_into   => {
               2 => 'store_peaks',
            },
        },
        {   -logic_name        => 'store_peaks',
            -module            => 'Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::StorePeaks',
            -analysis_capacity => 10,
        },

        {   -logic_name  => 'done_call_peaks',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
    ];
}

1;

