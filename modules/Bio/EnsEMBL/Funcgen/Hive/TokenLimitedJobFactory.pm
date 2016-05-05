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
  
  my $input_job = $self->input_job;
  
  $input_job->autoflow(0);
  
  my $this_logic_name  = $input_job->analysis->logic_name;
  my $this_analysis_id = $input_job->analysis->dbID;
  
  print "This logic name:  " . $this_logic_name . "\n";
  print "This analysis_id: " . $this_analysis_id . "\n";
  
  my $sql = qq(select from_analysis_id from dataflow_rule join dataflow_target on (dataflow_rule_id=source_dataflow_rule_id) where to_analysis_url = ? and dataflow_rule.from_analysis_id != ?);
  
  my $sth = $self->dbc->prepare($sql);
  
  $sth->bind_param(1, $this_logic_name);
  $sth->bind_param(2, $this_analysis_id);
  
  $sth->execute;
  
  my $x = $sth->fetchrow_arrayref();
  
  die() unless(ref $x eq 'ARRAY');
  die("Couldn't find job pool!") unless(@$x==1);
  
  my $pool_analysis_id = $x->[0];
  
  print "The analysis_id of the job pool is: " . $pool_analysis_id . "\n";
  
  my $jobs_from_pool = $input_job->adaptor->fetch_all_by_analysis_id_status(
    $pool_analysis_id, 
    'READY'
  );
  
  # Pool is empty
  if (!@$jobs_from_pool) {
    print "There are no available job in the job pool.\n";
    return;
  }
  
  print "Found " . @$jobs_from_pool . " jobs with status 'READY' in the job pool.\n";

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



