package Bio::EnsEMBL::Funcgen::PipeConfig::ChIPSeqAnalysis::ConvertSignalToBed;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Funcgen::PipeConfig::ChIPSeqAnalysis::Base';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub pipeline_analyses {
    my $self = shift;
    return [
       {   -logic_name  => 'start_convert_signal_to_bed',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               'MAIN->A' => 'convert_signal_to_bed',
               'A->MAIN' => 'done_convert_signal_to_bed',
            },
        },
        {   -logic_name  => 'convert_signal_to_bed',
            -parameters  => {
              convert_controls => 0,
            },
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::ConvertBamToBed',
        },
        {   -logic_name  => 'done_convert_signal_to_bed',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
    ];
}

1;

