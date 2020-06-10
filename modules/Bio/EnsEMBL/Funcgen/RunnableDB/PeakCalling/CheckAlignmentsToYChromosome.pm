package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::CheckAlignmentsToYChromosome;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;

use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw ( :all );

sub run {

  my $self = shift;
  
  my $species       = $self->param_required('species');
  my $plan          = $self->param_required('execution_plan');
  my $data_root_dir = $self->param_required('data_root_dir');

  #This analysis will only run on Human
  if ($species ne 'homo_sapiens') {
    return;
  }

  use Bio::EnsEMBL::Funcgen::PeakCallingPlan::ExecutionPlanUtils qw (
        lock_execution_plan
        resolve_nonterminal_symbols
  );
  my $plan_expanded = resolve_nonterminal_symbols($plan);
  lock_execution_plan($plan_expanded);
  lock_execution_plan($plan);
  
  print Dumper($plan_expanded);

  my $aligned_to_gender = $plan_expanded->{'input'}->{'to_gender'};
  my $bam_file          = $plan_expanded->{'input'}->{'output'}->{'real'};
  
  my $check_bam_file = $data_root_dir . '/' . $bam_file;
  
  my $alignments_to_y_chromosome_expected = $aligned_to_gender eq 'male';

  if ($alignments_to_y_chromosome_expected) {
    $self->say_with_header("ok, alignments to Y chromosome are not an error.", 1);
    return;
  }
  
  my $cmd = "samtools view $check_bam_file | cut -f 3 | uniq | grep Y | head";
  $self->say_with_header("Running: " . $cmd, 1);
  
  my $output = `$cmd`;
  
  if ($output) {
    die("Found alignments to Y chromosome!");
  }
  
  $self->say_with_header("ok, no alignments to Y chromosome found.", 1);
  
  return;
}

1;
