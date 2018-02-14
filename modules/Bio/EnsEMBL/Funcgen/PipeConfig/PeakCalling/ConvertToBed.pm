package Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::ConvertToBed;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::Base';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub pipeline_analyses {
    my $self = shift;
    return [
       {   -logic_name  => 'start_convert_to_bed',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               'MAIN->A' => 'convert_to_bed',
               'A->MAIN' => 'done_convert_to_bed',
            },
        },
        {   -logic_name  => 'convert_to_bed',
            -parameters  => {
              convert_controls => 1,
            },
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::ConvertBamToBed',
        },
        {   -logic_name  => 'done_convert_to_bed',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
    ];
}

1;

