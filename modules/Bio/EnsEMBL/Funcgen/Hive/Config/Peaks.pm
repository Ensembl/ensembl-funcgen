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
      -logic_name => 'PreprocessAlignments',
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::CollectionWriter',
      -flow_into => {
	2   => [ 
	  WHEN(
	    '#feature_set_analysis_logic_name# eq "SWEmbl_R015"'      => { 'CallDefaultPeaks' => INPUT_PLUS() },
	    '#feature_set_analysis_logic_name# eq "ccat_histone"'     => { 'CallBroadPeaks'   => INPUT_PLUS() },
	    '#feature_set_analysis_logic_name# eq "SWEmbl_R0025"'     => { 'CallTightPeaks'   => INPUT_PLUS() },
	    '#feature_set_analysis_logic_name# eq "SWEmbl_R0005_IDR"' => { 'CallIDRPeaks'     => INPUT_PLUS() },
	  ),
	  'WriteBigWig'
	],
        7 => { 
            ':////accu?file_to_delete=[]' => { 
              'file_to_delete' => '#file_to_delete#'
            } 
          }
      },
    },
    {
     -logic_name    => 'WriteBigWig',
     -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunWiggleTools',
     -parameters    => { mode => 'RPKM' },
     -rc_name => 'normal_30GB_2cpu',
    },
    {
      # Normal peaks
      -logic_name    => 'CallDefaultPeaks',  
      -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunPeaks',
      -rc_name => 'normal_5GB_2cpu_monitored',
    },
    {
      -logic_name    => 'CallBroadPeaks',
      -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunPeaks', 
      -parameters    => {
	process_file_types => ['significant_region'],
	CCAT_parameters => {
	  -chr_file => $self->o('chromosome_file')
	}
      },
      -rc_name => 'normal_monitored_4GB',
    },
    {
      # Tight peaks 
      -logic_name    => 'CallTightPeaks',
      -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunPeaks', 
      -rc_name => 'normal_5GB_2cpu_monitored',
    },
    {
      # Permissive peaks filtered by IDR max peaks value
      -logic_name    => 'CallIDRPeaks',
      -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunPeaks', 
      -parameters    => {
	filter_max_peaks => 1
      },
      -rc_name => 'normal_5GB_2cpu_monitored',
    },
   ];
}

1;
