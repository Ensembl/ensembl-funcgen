package Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::SplitByIdrStrategy;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::IDRStrategy;

use constant {
  BRANCH_SKIP_IDR                         => 2,
  BRANCH_RUN_IDR_ON_TECHNICAL_REPLICATES  => 3,
  BRANCH_RUN_IDR_ON_BIOLOGICAL_REPLICATES => 4,
};

my $branch_map = {

  Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::IDRStrategy
    ->SKIP_IDR => BRANCH_SKIP_IDR,

  Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::IDRStrategy
    ->RUN_IDR_ON_TECHNICAL_REPLICATES => BRANCH_RUN_IDR_ON_TECHNICAL_REPLICATES,

  Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::IDRStrategy
    ->RUN_IDR_ON_BIOLOGICAL_REPLICATES => BRANCH_RUN_IDR_ON_BIOLOGICAL_REPLICATES,

};

sub run {

  my $self           = shift;
  my $species        = $self->param_required('species');
  my $execution_plan = $self->param_required('execution_plan');
  
  my $idr_strategy = $execution_plan
    ->{call_peaks}
    ->{run_idr}
    ->{type};
  
  if (! defined $idr_strategy) {
    die(
      "Can't find idr strategy in execution plan!\n" 
      . Dumper($execution_plan)
    );
  }
  
  my $branch = $branch_map->{$idr_strategy};
  
  if (! defined $branch) {
    die("Unknown idr strategy $idr_strategy!");
  }
  $self->dataflow_output_id( 
    {
      'execution_plan' => $execution_plan,
      'species'        => $species,
    }, 
    $branch
  );
  return;
}
1;
