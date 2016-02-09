package Bio::EnsEMBL::Funcgen::Hive::Config::BaseDB;

use strict;
use warnings;
use Data::Dumper;
use base qw(Bio::EnsEMBL::Funcgen::Hive::Config::Base);

sub default_options {
  my $self = $_[0];  
  
  return {
    %{$self->SUPER::default_options},
	dnadb_pass          => $self->o('ENV', 'DNADB_PASS'),
	pass                => $self->o('ENV', 'DB_PASS'),
	dnadb_port          => undef,
	port                => undef,

	# These can probably go:
	#
	ssh                 => undef, #Connect to DBs using ssh(use in Importer)
	result_set_only    => 0, #why is this 0 rather than undef?
   };
}

sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
      %{$self->SUPER::pipeline_wide_parameters},
      dnadb   => {
         -dnadb_host   => $self->o('dnadb_host'),
         -dnadb_pass   => $self->o('dnadb_pass'),
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
    };
}

1;
