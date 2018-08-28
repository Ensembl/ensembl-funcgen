package Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::ConvertSignalToBed;

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
               'MAIN->A' => 'convert_to_bed_multiplex',
               'A->MAIN' => 'done_convert_signal_to_bed',
            },
        },
       {   -logic_name  => 'convert_to_bed_multiplex',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               MAIN => [
                'convert_signal_to_bed',
                'convert_control_to_bed',
               ]
            },
        },
        {   -logic_name  => 'convert_signal_to_bed',
            -parameters  => {
              convert_controls => 0,
            },
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::ConvertBamToBed',
#             -flow_into   => {
#                MEMLIMIT => 'convert_signal_to_bed_himem',
#             },
        },
        {   -logic_name  => 'convert_control_to_bed',
            -parameters  => {
              convert_controls => 1,
            },
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::ConvertBamToBed',
        },
#         {   -logic_name  => 'convert_signal_to_bed_himem',
#             -parameters  => {
#               convert_controls => 0,
#             },
#             -rc_name     => '8Gb_job',
#             -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::ConvertBamToBed',
#         },
        {   -logic_name  => 'done_convert_signal_to_bed',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
    ];
}

1;

