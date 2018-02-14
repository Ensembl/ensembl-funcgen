package Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::CallPeaks;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::Base';
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
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::CallPeaks',
            -parameters => {
              tempdir => '#tempdir_peak_calling#/#species#/peak_calling'
            },
            -flow_into   => {
               2        => 'store_peaks',
               MEMLIMIT => 'call_peaks_himem',
            },
        },
        {   -logic_name  => 'call_peaks_himem',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::CallPeaks',
            -parameters => {
              tempdir => '#tempdir_peak_calling#/#species#/peak_calling'
            },
            -rc_name     => '8Gb_job',
            -flow_into   => {
               2 => 'store_peaks',
            },
        },
        {   -logic_name        => 'store_peaks',
            -module            => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::StorePeaks',
            -analysis_capacity => 10,
        },
        {   -logic_name  => 'done_call_peaks',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               1 => 'fix_known_problems',
            },
        },
        {   -logic_name => 'fix_known_problems',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -analysis_capacity => 1,
            -parameters => {
            
                peak_calling => '#expr( #execution_plan#->{call_peaks}->{name} )expr#',
                tempdir      => '#tempdir_peak_calling#/#species#/peak_calling/#peak_calling#',
                
                cmd => qq( trim_peaks_to_seq_region_boundaries.pl )
                . qq( --species      #species#      )
                . qq( --registry     #reg_conf#     )
                . qq( --peak_calling #peak_calling# )
                . qq( --tempdir      #tempdir#      )
            },
        },
    ];
}

1;

