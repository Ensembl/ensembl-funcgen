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
            -analysis_capacity => 100,
            -rc_name     => '4Gb_job_2h',
            -flow_into   => {
               2        => {
                'store_peaks' => undef,
               },
               MEMLIMIT => {
                'call_peaks_himem' => undef,
               },
               RUNLIMIT => {
                'store_peaks' => INPUT_PLUS({
                      peak_calling_succeeded => 0,
                      error_message          => 'Job exceeded maximal run time',
                  })
                },
            },
        },
        {   -logic_name  => 'call_peaks_himem',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::CallPeaks',
            -parameters => {
              tempdir => '#tempdir_peak_calling#/#species#/peak_calling'
            },
            -analysis_capacity => 50,
            -rc_name     => '32Gb_job_2h',
            -flow_into   => {
               2 => { 
                  'store_peaks' => undef,
               },
               RUNLIMIT => {
                'store_peaks' => INPUT_PLUS({
                      peak_calling_succeeded => 0,
                      error_message          => 'Job exceeded maximal run time',
                  })
                },
            },
        },
        {   -logic_name        => 'store_peaks',
            -module            => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::StorePeaks',
            -rc_name           => '16Gb_job',
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

