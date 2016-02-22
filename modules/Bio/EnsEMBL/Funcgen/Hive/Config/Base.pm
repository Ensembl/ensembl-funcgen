package Bio::EnsEMBL::Funcgen::Hive::Config::Base;

use strict;
use warnings;
use base qw(Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf);

sub default_options {
  my $self = shift;
  
  return {
      %{$self->SUPER::default_options},
      dnadb_pass        => $self->o('ENV', 'DNADB_PASS'),
      pass              => $self->o('ENV', 'DB_PASS'),
      use_tracking_db   => 1,
      work_root_dir     => $self->o('data_root_dir').'/output/'.$self->o('pipeline_name'),
      hive_output_dir   => $self->o('data_root_dir').'/output/'.$self->o('pipeline_name').'/hive_debug',
      port              => 3306,
      dnadb_port        => 3306,
      pipeline_name     => 'ersa',
   };
}

sub pipeline_wide_parameters {
  my ($self) = @_;

  my $dnadb_pass = $self->o('dnadb_pass');
  # This is code for no password on the command line.
  #
  if (($dnadb_pass eq "''") || ($dnadb_pass eq "\"\"")) {
    $dnadb_pass = undef;
  }
  
  return {
    %{$self->SUPER::pipeline_wide_parameters},

    dnadb   => {
	-dnadb_host   => $self->o('dnadb_host'),
	-dnadb_pass   => $dnadb_pass,
	-dnadb_port   => $self->o('dnadb_port'),
	-dnadb_user   => $self->o('dnadb_user'),
	-dnadb_name   => $self->o('dnadb_name'),
    },
    out_db  => {
	-host   => $self->o('host'),
	-port   => $self->o('port'),
	-user   => $self->o('user'),
	-pass   => $self->o('pass'),
	-dbname => $self->o('dbname'),
    },
    pipeline_name => $self->o('pipeline_name'),

    species        => $self->o('species'),
    assembly       => $self->o('assembly'),
    data_root_dir     => $self->o('data_root_dir'),
    work_root_dir     => $self->o('work_root_dir'),
    hive_output_dir => $self->o('hive_output_dir'),
    use_tracking_db => $self->o('use_tracking_db'),

    default_gender => 'male',
  };
}

sub resource_classes {
  my $self = shift;
  return {
     default                 => { 'LSF' => '' },    
     normal_2GB              => { 'LSF' => ' -M2000 -R"select[mem>2000] rusage[mem=2000] span[hosts=1]"' },
     normal_monitored        => { 'LSF' => "" },
     normal_high_mem         => { 'LSF' => ' -M5000 -R"select[mem>5000] rusage[mem=5000] span[hosts=1]"' },
     normal_high_mem_2cpu    => { 'LSF' => ' -n2 -M5000 -R"select[mem>5000] rusage[mem=5000] span[hosts=1]"' },
     normal_monitored_2GB    => {'LSF' => " -M2000 -R\"select[mem>2000]".
                                                " rusage[mem=2000] span[hosts=1]\"" },
     normal_monitored_4GB    => {'LSF' => " -M4000 -R\"select[mem>4000] rusage[mem=4000] span[hosts=1]\"" },
     normal_4GB_2cpu         => {'LSF' => " -n2 -M4000 -R\"select[mem>4000] rusage[mem=4000] span[hosts=1]\"" },
     normal_monitored_8GB    => {'LSF' => " -M8000 -R\"select[mem>8000] rusage[mem=8000]\"" },
     normal_monitored_8GB_2cpu => {'LSF' => " -n2 -M8000 -R\"select[mem>8000] rusage[mem=8000]\"" },   
     normal_monitored_16GB   => {'LSF' => " -M16000 -R\"select[mem>16000] rusage[mem=16000]\"" }, 
     normal_16GB_2cpu        => {'LSF' => ' -n2 -M16000 -R"select[mem>16000] rusage[mem=16000] span[hosts=1]"' },
     normal_20GB_2cpu        => {'LSF' => ' -n2 -M20000 -R"select[mem>20000] rusage[mem=20000] span[hosts=1]"' }, 
     normal_25GB_2cpu        => {'LSF' => ' -n2 -M25000 -R"select[mem>25000] rusage[mem=25000] span[hosts=1]"' }, 
     normal_30GB_2cpu        => {'LSF' => ' -n2 -M30000 -R"select[mem>30000] rusage[mem=30000] span[hosts=1]"' },
     normal_30GB_3cpu        => {'LSF' => ' -n3 -M30000 -R"select[mem>30000] rusage[mem=30000] span[hosts=1]"' },
     '64GB_3cpu'             => {'LSF' => ' -n3 -M64000 -R"select[mem>64000] rusage[mem=64000] span[hosts=1]"' },
     normal_10gb_monitored   => {'LSF' => " -M10000 -R\"select[mem>10000] rusage[mem=10000] span[hosts=1]\"" },
     normal_5GB_2cpu_monitored => {'LSF' => " -n2 -M5000 -R\"select[mem>5000] rusage[mem=5000] span[hosts=1]\"" },
     normal_10gb             => { 'LSF' => ' -M10000 -R"select[mem>10000] rusage[mem=10000] span[hosts=1]"' },
     long_monitored          => { 'LSF' => "-q long " },
     long_high_mem           => { 'LSF' => '-q long -M4000 -R"select[mem>4000] rusage[mem=4000] span[hosts=1]"' },
     long_monitored_high_mem => { 'LSF' => "-q long -M4000 -R\"select[mem>4000] rusage[mem=4000] span[hosts=1]\"" },
    };
}

1;
