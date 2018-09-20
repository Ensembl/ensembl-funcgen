package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::RunPermissiveSWEmbl;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;

use constant {
  BRANCH_OUTPUT => 2,
};

sub run {

  my $self = shift;
  
  my $data_root_dir = $self->param_required('data_root_dir');
  my $tempdir       = $self->param_required('tempdir');
  my $plan          = $self->param_required('execution_plan');
  
  use Bio::EnsEMBL::Funcgen::PeakCallingPlan::ExecutionPlanUtils qw (
        lock_execution_plan
        resolve_nonterminal_symbols
  );
  my $plan_expanded = resolve_nonterminal_symbols($plan);
  lock_execution_plan($plan_expanded);

  my $bam_file_no_duplicates = $plan_expanded
    ->{output}
    ->{real}
  ;
  my $alignment_name = $plan_expanded->{name};

  my $job_output = $tempdir . '/' . $bam_file_no_duplicates . '.swembl.txt';
  
  use File::Basename qw( dirname basename );
  my $dirname = dirname($job_output);

  use File::Path qw(make_path remove_tree);
  make_path($dirname);

  my $swembl_parameters = "-f 150 -m 8 -p 0.04 -P 0.5 -d 70 ";
  my $cmd = "SWEmbl -F $swembl_parameters -i $data_root_dir/$bam_file_no_duplicates -o $job_output";
  my $error_occurred = $self->run_system_command($cmd);
  
  if ($error_occurred) {
    $self->throw("There was an error running:\n$cmd");
  }
  $self->dataflow_output_id( 
    {
      'permissive_peak_calling' => {
          'peak_file'      => $job_output,
          'alignment_name' => $alignment_name,
      }
    }, 
    BRANCH_OUTPUT
  );

}

1;
