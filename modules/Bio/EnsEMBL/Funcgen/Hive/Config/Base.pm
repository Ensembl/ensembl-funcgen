=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 NAME

    Bio::EnsEMBL::Funcgen::Hive::Config::Base;

=head1 CONTACT

    Please contact http://lists.ensembl.org/mailman/listinfo/dev mailing list with questions/suggestions.

=cut
package Bio::EnsEMBL::Funcgen::Hive::Config::Base;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf);

sub default_options {
  my $self = shift;
  
  return {
      %{$self->SUPER::default_options},
      pipeline_name => 'ersa4ever',
   };
}

sub beekeeper_extra_cmdline_options {
    my ($self) = @_;
    return '-reg_conf ' . $self->o('reg_conf') . ' -keep_alive -can_respecialize 1 -sleep 0.1';
}

sub pipeline_wide_parameters {
  my ($self) = @_;
  
  my $species  = $self->o('species');
  my $reg_conf = $self->o('reg_conf');
  
  my $second_pass = $species!~/^#:subst/ && $reg_conf!~ /^#:subst/;
  
  my $details_from_registry = {};
  
  if ($second_pass) {
    my $reg_conf = $self->o('reg_conf');
    my $species  = $self->o('species');
    $details_from_registry = fetch_details_from_registry($species, $reg_conf);
  }
  
  return {
    %{$self->SUPER::pipeline_wide_parameters},

    dnadb            => $details_from_registry->{dnadb},
    out_db           => $details_from_registry->{outdb},
    pipeline_name    => $self->o('pipeline_name'),
    ensembl_release_version => $self->o('ensembl_release_version'),

    species          => $details_from_registry->{species},
    assembly         => $details_from_registry->{assembly},
    data_root_dir    => $self->o('data_root_dir'),
    work_root_dir    => $self->o('tempdir') . '/chip_seq_analysis',
    use_tracking_db  => 1,
  };
}

sub fetch_details_from_registry {

  my $species  = shift;
  my $registry = shift;
  
  use Bio::EnsEMBL::Registry;
  Bio::EnsEMBL::Registry->load_all($registry);

  my $funcgen_db_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
  my $core_db_adaptor    = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'core');
  
  my $coordsystem_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'coordsystem');

  my $default_chromosome_coordsystem = $coordsystem_adaptor->fetch_by_name('chromosome');
  my $default_assembly = $default_chromosome_coordsystem->version;
  
  my $dnadb_host = $core_db_adaptor->dbc->host;
  my $dnadb_pass = $core_db_adaptor->dbc->pass;
  my $dnadb_port = $core_db_adaptor->dbc->port;
  my $dnadb_user = $core_db_adaptor->dbc->user;
  my $dnadb_name = $core_db_adaptor->dbc->dbname;
  
  my $host   = $funcgen_db_adaptor->dbc->host;
  my $pass   = $funcgen_db_adaptor->dbc->pass;
  my $port   = $funcgen_db_adaptor->dbc->port;
  my $user   = $funcgen_db_adaptor->dbc->user;
  my $dbname = $funcgen_db_adaptor->dbc->dbname;

  my $dnadb = {
    "-dnadb_host" => "${dnadb_host}",
    "-dnadb_name" => "$dnadb_name",
    "-dnadb_pass" => "$dnadb_pass",
    "-dnadb_port" => $dnadb_port,
    "-dnadb_user" => "$dnadb_user"
  };
  my $outdb = {
    "-dbname" => "$dbname",
    "-driver" => "mysql",
    "-host" => "$host",
    "-pass" => "$pass",
    "-port" => $port,
    "-user" => "$user"
  };

  use Hash::Util qw( lock_hash );
  my $registry_details = {

    dnadb_host   => $dnadb_host,
    dnadb_pass   => $dnadb_pass,
    dnadb_port   => $dnadb_port,
    dnadb_user   => $dnadb_user,
    dnadb_name   => $dnadb_name,
    
    species => $core_db_adaptor->species,
    
    dnadb => $dnadb,
    outdb => $outdb,

    host   => $host,
    pass   => $pass,
    port   => $port,
    user   => $user,
    name   => $dbname,

    assembly => $default_assembly
  };
  
  return $registry_details;
}

sub resource_classes {
  my $self = shift;
  return {
     default                 => { 'LSF' => '' },    
     normal_2GB              => { 'LSF' => ' -q production-rh7 -M2000 -R"select[mem>2000] rusage[mem=2000] span[hosts=1]"' },
     normal_monitored        => { 'LSF' => "" },
     normal_high_mem         => { 'LSF' => ' -q production-rh7 -M5000 -R"select[mem>5000] rusage[mem=5000] span[hosts=1]"' },
     normal_high_mem_2cpu    => { 'LSF' => ' -q production-rh7 -n2 -M5000 -R"select[mem>5000] rusage[mem=5000] span[hosts=1]"' },
     normal_monitored_2GB    => {'LSF' => " -q production-rh7 -M2000 -R\"select[mem>2000]".
                                                " rusage[mem=2000] span[hosts=1]\"" },
     normal_monitored_4GB    => {'LSF' => " -q production-rh7 -M4000 -R\"select[mem>4000] rusage[mem=4000] span[hosts=1]\"" },
     normal_4GB_2cpu         => {'LSF' => " -q production-rh7 -n2 -M4000 -R\"select[mem>4000] rusage[mem=4000] span[hosts=1]\"" },
     normal_monitored_8GB    => {'LSF' => " -q production-rh7 -M8000 -R\"select[mem>8000] rusage[mem=8000]\"" },
     normal_monitored_8GB_2cpu => {'LSF' => " -q production-rh7 -n2 -M8000 -R\"select[mem>8000] rusage[mem=8000]\"" },   
     normal_monitored_16GB   => {'LSF' => " -q production-rh7 -M16000 -R\"select[mem>16000] rusage[mem=16000]\"" }, 
     normal_16GB_2cpu        => {'LSF' => ' -q production-rh7 -n2 -M16000 -R"select[mem>16000] rusage[mem=16000] span[hosts=1]"' },
     normal_20GB_2cpu        => {'LSF' => ' -q production-rh7 -n2 -M20000 -R"select[mem>20000] rusage[mem=20000] span[hosts=1]"' }, 
     normal_25GB_2cpu        => {'LSF' => ' -q production-rh7 -n2 -M25000 -R"select[mem>25000] rusage[mem=25000] span[hosts=1]"' }, 
     normal_30GB_2cpu        => {'LSF' => ' -q production-rh7 -n2 -M30000 -R"select[mem>30000] rusage[mem=30000] span[hosts=1]"' },
     normal_30GB_3cpu        => {'LSF' => ' -q production-rh7 -n3 -M30000 -R"select[mem>30000] rusage[mem=30000] span[hosts=1]"' },
     '64GB_3cpu'             => {'LSF' => ' -q production-rh7 -n3 -M64000 -R"select[mem>64000] rusage[mem=64000] span[hosts=1]"' },
     #'10gb_1cpu_staggered'   => {'LSF' => q(-E 'sleep $(echo "$LSB_JOBINDEX * 1" | bc)' -M10000 -R"select[mem>10000] rusage[mem=10000] span[hosts=1]") },
     '10gb_1cpu'             => {'LSF' => q( -q production-rh7     -M10000 -R"select[mem>10000] rusage[mem=10000] span[hosts=1]") },
     '10gb_2cpu'             => {'LSF' => q( -q production-rh7 -n2 -M10000 -R"select[mem>10000] rusage[mem=10000] span[hosts=1]") },
      
     normal_5GB_2cpu_monitored => {'LSF' => " -q production-rh7 -n2 -M5000 -R\"select[mem>5000] rusage[mem=5000] span[hosts=1]\"" },
     normal_10gb             => { 'LSF' => ' -q production-rh7 -M10000 -R"select[mem>10000] rusage[mem=10000] span[hosts=1]"' },
     long_monitored          => { 'LSF' => " -q production-rh7 " },
     long_high_mem           => { 'LSF' => ' -q production-rh7 -M4000 -R"select[mem>4000] rusage[mem=4000] span[hosts=1]"' },
     long_monitored_high_mem => { 'LSF' => " -q production-rh7 -M4000 -R\"select[mem>4000] rusage[mem=4000] span[hosts=1]\"" },
    };
}

1;
