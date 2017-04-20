package Bio::EnsEMBL::Funcgen::Hive::Config::Peaks;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('Bio::EnsEMBL::Funcgen::Hive::Config::Base');

sub pipeline_analyses {
  my $self = shift;

  return 
   [
    {
      -logic_name => 'index_bam_files',
      -flow_into => {
	2   => [ 
	  'start_peak_calling',
	  'write_bigwig'
	],
        7 => { 
            ':////accu?file_to_delete=[]' => { 
              'file_to_delete' => '#file_to_delete#'
            } 
          }
      },
    },
    {
     -logic_name => 'write_bigwig',
     -module     => 'Bio::EnsEMBL::Funcgen::Hive::RunWiggleTools',
     -rc_name    => 'normal_30GB_2cpu',
    },
    {
      -logic_name => 'start_peak_calling',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      -flow_into => {
        'MAIN->A' => [
          WHEN(
              '#logic_name# eq "SWEmbl_R015"'      => 'call_default_peaks',
              '#logic_name# eq "ccat_histone"'     => 'call_broad_peaks_start',
              '#logic_name# eq "SWEmbl_R0025"'     => 'call_tight_peaks',
              '#logic_name# eq "SWEmbl_R0005_IDR"' => 'call_idr_peaks',
            ),
        ],
        'A->MAIN' => [ 'done_peak_calling' ]
      }
    },
    {
      -logic_name => 'done_peak_calling',
      -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
    },
    {
      # Normal peaks
      -logic_name    => 'call_default_peaks',  
      -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunPeaks',
      -rc_name => 'normal_5GB_2cpu_monitored',
    },

    {   -logic_name => 'call_broad_peaks_start',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        -flow_into => {
          'MAIN->A' => 'convert_bam_to_bed',
          'A->MAIN' => 'call_broad_peaks',
        },
    },
    {
      -logic_name    => 'convert_bam_to_bed',
      -module        => 'Bio::EnsEMBL::Funcgen::Hive::ConvertBamToBed', 
    },
    {
      -logic_name    => 'call_broad_peaks',
      -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunPeaks', 
      -parameters    => {
        process_file_types => ['significant_region'],
      },
      -rc_name => 'normal_monitored_4GB',
    },

    {
      # Tight peaks 
      -logic_name    => 'call_tight_peaks',
      -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunPeaks', 
      -rc_name => 'normal_5GB_2cpu_monitored',
    },
    {
      # Permissive peaks filtered by IDR max peaks value
      -logic_name    => 'call_idr_peaks',
      -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunPeaks', 
      -parameters    => {
	filter_max_peaks => 1
      },
      -rc_name => 'normal_5GB_2cpu_monitored',
    },
   ];
}

1;
