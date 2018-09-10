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
      {   -logic_name  => 'start_populate_read_file_stats',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -flow_into => {
              'MAIN->A' => 'find_read_file_entries',
              'A->MAIN' => 'generate_read_file_report',
          },
      },
      {   -logic_name  => 'find_read_file_entries',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
          -parameters => {
              db_conn    => 'funcgen:#species#',
              inputquery => '
                select 
                  read_file_id 
                from 
                  read_file 
                where 
                  file_size = 0 
                  or read_length = 0 
                  or file_size is null 
                  or read_length is null
              ',
          },
          -flow_into => {
              2 => { 'populate_read_file_stats',  INPUT_PLUS },
          },
      },
      {   -logic_name => 'populate_read_file_stats',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          #-analysis_capacity => 50,
          -parameters => {
            cmd => qq( compute_read_file_statistics.pl )
              . qq( --species         #species#        )
              . qq( --registry        #reg_conf#       )
              . qq( --read_file_id    #read_file_id#   )
          },
      },
      {   -logic_name => 'generate_read_file_report',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 
                qq( generate_read_file_report.pl                 )
              . qq(   --species          #species#               )
              . qq(   --registry         #reg_conf#              )
              . qq(   --output_directory #reports_dir#/#species# )
          },
      },
    ]
}

1;
