package Bio::EnsEMBL::Funcgen::Hive::Config::Collections;

use strict;
use warnings;
use base ('Bio::EnsEMBL::Funcgen::Hive::Config::BaseSequenceAnalysis');

sub pipeline_analyses {
  my $self = shift;

  return [
     {
      -logic_name => 'PreprocessAlignments',
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::CollectionWriter',
      -flow_into => { 2 => ['WriteBigWig'] }, 
      -analysis_capacity => 50,
      -rc_name => 'normal_2GB',
    },
    {
     -logic_name    => 'WriteBigWig',
     -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunWiggleTools',
     -parameters    => { mode => 'RPKM' },
     -analysis_capacity => 100,
     -rc_name => 'normal_30GB_2cpu',
    },
  ];
}

1;
