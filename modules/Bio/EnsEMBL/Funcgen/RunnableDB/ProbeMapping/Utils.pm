package Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::Utils;

use warnings;
use strict;


use base qw( Exporter );
use vars qw( @EXPORT_OK );

@EXPORT_OK = qw(
  create_probe_adaptor
  create_funcgen_adaptor
  create_dna_db_adaptor
  create_db_url_from_dbc
  create_db_server_url_from_dbc
  create_db_url_from_dba_hash
  create_db_server_url_from_dba_hash
  create_dna_db_params_from_funcgen_hash
);

sub create_db_url_from_dba_hash {

  my $dba_hash = shift;
  
  # Making a detour via DBConnection. DBConnection set defaults,
  # if they are missing from the hash, that way we get Ensembl typical 
  # behaviour for missing values.
  #
  my $dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(%$dba_hash);
  
  my $url = create_db_url_from_dbc($dbc);
  
  my @query_parameters;
  if (exists $dba_hash->{-group}) {
    push @query_parameters, 'group=' . $dba_hash->{-group};
  }
  if (exists $dba_hash->{-species}) {
    push @query_parameters, 'species=' . $dba_hash->{-species};
  }
  $url .= '?' . join '&', @query_parameters;
  
  return $url;
}

sub create_db_server_url_from_dba_hash {

  my $dba_hash = shift;
  my $dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(%$dba_hash);
  return create_db_server_url_from_dbc($dbc);
}

sub create_db_url_from_dbc {
  my $dbc = shift;
  return create_db_server_url_from_dbc($dbc) . $dbc->dbname;

}
sub create_db_server_url_from_dbc {
  my $dbc = shift;
  return $dbc->driver . '://'.$dbc->user.':'.$dbc->password.'@'.$dbc->host.':'.$dbc->port . '/';
}

sub create_probe_adaptor {
  my $funcgen_dba = create_funcgen_adaptor(@_);  
  return $funcgen_dba->get_ProbeAdaptor;
}

sub create_dna_db_adaptor {
  my $dba_hash = shift;
  
  use Bio::EnsEMBL::DBSQL::DBAdaptor;
  my $dnadb =  Bio::EnsEMBL::DBSQL::DBAdaptor->new(
    create_dna_db_params_from_funcgen_hash($dba_hash)
  );
  return $dnadb;
}

sub create_dna_db_params_from_funcgen_hash {
  my $funcgen_hash = shift;
  my @params = (
      -dbname          => $funcgen_hash->{'-dnadb_name'},
      -host            => $funcgen_hash->{'-dnadb_host'},
      -port            => $funcgen_hash->{'-dnadb_port'},
      -user            => $funcgen_hash->{'-dnadb_user'},
      -species         => $funcgen_hash->{'-species'}, 
  );
  return @params;
}

sub create_funcgen_adaptor {
  my $dba_hash = shift;
  
  use Bio::EnsEMBL::DBSQL::DBAdaptor;
  my $dnadb = create_dna_db_adaptor($dba_hash);

  delete $dba_hash->{'-dnadb_name'};
  delete $dba_hash->{'-dnadb_host'};
  delete $dba_hash->{'-dnadb_port'};
  delete $dba_hash->{'-dnadb_user'};

  use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
  my $funcgen_dbadaptor = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
    (
      %$dba_hash, 
      -dnadb => $dnadb,
    );
  return $funcgen_dbadaptor;
}

1;

