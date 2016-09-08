package Bio::EnsEMBL::Funcgen::Hive::Ftp::DbconnForSpecies;

use strict;
use Data::Dumper;
use base ('Bio::EnsEMBL::Hive::Process');

sub run {
  my $self    = shift;
  my $species = $self->param('species');
  my $group   = $self->param('group');

  use Bio::EnsEMBL::Registry;
  my $adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, $group );  
  my $dbc = Bio::EnsEMBL::Hive::DBSQL::DBConnection->new(-dbconn => $adaptor->dbc);
  
  $self->dataflow_output_id({
    url => $dbc->url
  }, 2);
}

1;
