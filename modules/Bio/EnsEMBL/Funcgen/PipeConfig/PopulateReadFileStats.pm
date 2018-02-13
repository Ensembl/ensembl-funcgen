package Bio::EnsEMBL::Funcgen::PipeConfig::PopulateReadFileStats;

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
      pipeline_name => 'populate_read_file_stats',
   };
}

sub pipeline_wide_parameters {
  my $self = shift;
  
  return {
    %{$self->SUPER::pipeline_wide_parameters},
    pipeline_name           => $self->o('pipeline_name'),
    reg_conf                => $self->o('reg_conf'),
  };
}

sub pipeline_analyses {
    my $self = shift;

    return [
      {   -logic_name  => 'start',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
          -parameters => {
              db_conn    => 'funcgen:#species#',
              inputquery => '
                select 
                  read_file_id
                from 
                  read_file
              ',
          },
          -flow_into => {
              2 => { 'populate_read_file_stats' => INPUT_PLUS },
          },
      },
      {   -logic_name => 'populate_read_file_stats',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -analysis_capacity => 50,
          -parameters => {
            cmd => qq( compute_read_file_statistics.pl )
              . qq( --species         #species#        )
              . qq( --registry        #reg_conf#       )
              . qq( --read_file_id    #read_file_id#   )
          },
      },
    ]
}

1;
