=head1 NAME

Bio::EnsEMBL::Funcgen::Hive::TokenLimitedJobFactory

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::Hive::TokenLimitedJobFactory;

use base ('Bio::EnsEMBL::Hive::Process');
use Data::Dumper;
use warnings;
use strict;

my $self_branch     = 1;
my $pipeline_branch = 2;

my $num_jobs_a_token_seeds = 1;

sub run {
  my $self = shift;
  
  my $dfr_adaptor = $self->db->get_DataflowRuleAdaptor;
  my $input_job = $self->input_job;
  
  $input_job->autoflow(0);
  
  my $this_logic_name  = $input_job->analysis->logic_name;
  my $this_analysis_id = $input_job->analysis->dbID;
  
  print "This logic name:  " . $this_logic_name . "\n";
  print "This analysis_id: " . $this_analysis_id . "\n";
  
  my $pool_analysis_id;

  my $data_flow_rule = $dfr_adaptor->fetch_all;
  RULE: foreach my $current_data_flow_rule (@$data_flow_rule) {
    
    my $job_pool_found = 
      $current_data_flow_rule->to_analysis_url eq $this_logic_name
      && $current_data_flow_rule->from_analysis_id != $this_analysis_id;
      
    if ($job_pool_found) {
      $pool_analysis_id = $current_data_flow_rule->from_analysis_id;
      last RULE;
    }
  }
  
  die("Couldn't find job pool!") unless($pool_analysis_id);
  
  my $jobs_from_pool = $input_job->adaptor->fetch_all_by_analysis_id_status(
    $pool_analysis_id, 
    'READY'
  );
  
  # Pool is empty
  if (!@$jobs_from_pool) {
    return;
  }
  my $count = $num_jobs_a_token_seeds;  
  JOB: foreach my $current_job (@$jobs_from_pool) {
  
    last JOB if ($count==0);
    $count--;
  
    my $input_id = $current_job->input_id;
    $self->dataflow_output_id($input_id, $pipeline_branch);
    $current_job->status('DONE');
    $current_job->adaptor->check_in_job($current_job);
  }
  
  $self->dataflow_output_id(
    { 
      make_me_unique => { time => '' . localtime(), pid => $$, host => $ENV{HOSTNAME} }
    }, 
    $self_branch
  );
  
  return;
}

1;



