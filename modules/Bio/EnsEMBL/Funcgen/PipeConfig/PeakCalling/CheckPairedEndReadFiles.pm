package Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::CheckPairedEndReadFiles;

use strict;
use warnings;
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base 'Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf';

sub beekeeper_extra_cmdline_options {
  my $self = shift;
  return '-reg_conf ' . $self->o('reg_conf') . ' -keep_alive -can_respecialize 1 -sleep 0.5';
}

sub default_options {
  my $self = shift;
  
  return {
      %{$self->SUPER::default_options},
      pipeline_name => 'check_paired_end_read_files',
   };
}

sub pipeline_wide_parameters {
  my $self = shift;
  
  return {
    %{$self->SUPER::pipeline_wide_parameters},
    pipeline_name => $self->o('pipeline_name'),
    reg_conf      => $self->o('reg_conf'),
  };
}

sub pipeline_analyses {
    my $self = shift;

    return [
      {   -logic_name  => 'start_check_paired_end_read_files',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -flow_into => {
              'MAIN->A' => 'find_execution_plans',
              'A->MAIN' => 'done_check_paired_end_read_files',
          },
      },
      {   -logic_name  => 'find_execution_plans',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
          -parameters => {
              db_conn    => 'funcgen:#species#',
              inputquery => 'select execution_plan_id from execution_plan',
          },
          -flow_into => {
              2 => { 'check_paired_end_read_files',  INPUT_PLUS },
          },
      },
      {   -logic_name => 'check_paired_end_read_files',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          #-analysis_capacity => 50,
          -parameters => {
            cmd => qq( check_paired_end_read_files.pl         )
              . qq( --species           #species#             )
              . qq( --registry          #reg_conf#            )
              . qq( --execution_plan_id #execution_plan_id#   )
          },
      },
      {   -logic_name => 'done_check_paired_end_read_files',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      },
    ]
}

1;
