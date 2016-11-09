package Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::ConnectionDetailsAsParameters;

use strict;
use base ('Bio::EnsEMBL::Hive::Process');

sub run {
    my $self = shift;

    my $species = $self->param('species');
    my $type    = $self->param('type');

    my $adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, $type);
    
    $self->dataflow_output_id( {
      dbname   => $adaptor->dbc->dbname,
      driver   => $adaptor->dbc->driver,
      host     => $adaptor->dbc->host,
      password => $adaptor->dbc->password,
      port     => $adaptor->dbc->port,
      username => $adaptor->dbc->username,
      species  => $adaptor->species,
    }, 1 );
}

1;





