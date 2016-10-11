package Bio::EnsEMBL::Funcgen::Hive::Ftp::FetchParametersFromDb;

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use base ('Bio::EnsEMBL::Hive::Process');

sub run {
  my $self = shift;

  my $species  = $self->param('species');

  my $coord_system_adaptor = Bio::EnsEMBL::Registry->get_adaptor( $species, 'Core', 'CoordSystem' );
  if (!$coord_system_adaptor) {
    die("Can't get coord system adaptor! Please configure your registry accordingly.")
  }
  my ($cs) = @{$coord_system_adaptor->fetch_all()};
  my $assembly = $cs->version();
  if (!$assembly) {
    die("Can't work out assembly for $species!")
  }

  $self->dataflow_output_id({
    species                    => $species,
    assembly                   => $assembly,
  }, 1);
}

1;
