package Bio::EnsEMBL::Funcgen::Hive::Ftp::JobFactoryAnnotatedFeaturesUsingIdLists;

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use base ('Bio::EnsEMBL::Hive::Process');

sub run {
  my $self = shift;

  my $temp_dir = $self->param('temp_dir');
  my $species  = $self->param('species');

  use File::Path qw(make_path remove_tree);
  make_path($temp_dir);

  my $sql = <<SQL
    select 
      distinct epigenome.production_name as epigenome_production_name, 
      feature_type.name as feature_type_name, 
      analysis.logic_name as analysis_logic_name 
    from 
      annotated_feature 
      join feature_set using (feature_set_id) 
      join epigenome using (epigenome_id) 
      join feature_type using (feature_type_id) 
      join analysis on (analysis.analysis_id=feature_set.analysis_id);
SQL
;

  my $coord_system_adaptor = Bio::EnsEMBL::Registry->get_adaptor( $species, 'Core', 'CoordSystem' );
  if (!$coord_system_adaptor) {
    die("Can't get coord system adaptor! Please configure your registry accordingly.")
  }
  my ($cs) = @{$coord_system_adaptor->fetch_all()};
  my $assembly = $cs->version();
  if (!$assembly) {
    die("Can't work out assembly for $species!")
  }
  
  my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'Funcgen' );
  my $dbc = $funcgen_adaptor->dbc;
  
  use Bio::EnsEMBL::Utils::SqlHelper;
  my $helper = Bio::EnsEMBL::Utils::SqlHelper->new(
    -DB_CONNECTION => $funcgen_adaptor->dbc
  );
  
  my @all_annotated_feature_batches;
  
  $helper->execute_no_return(
    -SQL      => $sql,
    -CALLBACK => sub {
    
      my $row = shift;

      my $epigenome_production_name = $row->[0];
      my $feature_type_name         = $row->[1];
      my $analysis_logic_name       = $row->[2];

      $self->dataflow_output_id({
        species                    => $species,
        assembly                   => $assembly,
        feature_type_name          => $feature_type_name,
        epigenome_production_name  => $epigenome_production_name,
        analysis_logic_name        => $analysis_logic_name,
      }, 2);
    }
  );
}

1;
