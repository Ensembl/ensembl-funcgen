package Bio::EnsEMBL::Funcgen::RunnableDB::Ftp::CreatePeakCallingDescriptionJobs;

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use base ('Bio::EnsEMBL::Hive::Process');

use constant {
  BRANCH_JOB_OUTPUT => 2,
};

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

  my $peak_calling_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'peakcalling');
  #my $peak_callings = $peak_calling_adaptor->fetch_all('LIMIT 20');
  my $peak_callings = $peak_calling_adaptor->fetch_all;
  
  my $epigenome_adaptor    = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'epigenome');
  my $feature_type_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'featuretype');
  my $analysis_adaptor     = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'analysis');
  
  # Find all unique pairs
  my %unique;
  foreach my $peak_calling (@$peak_callings) {
  
    my $epigenome_id    = $peak_calling->epigenome_id;
    my $feature_type_id = $peak_calling->feature_type_id;
    
    $unique{$epigenome_id}{$feature_type_id} = $peak_calling->analysis_id;
  }
  my @all_epigenome_ids = keys %unique;
  
  # Create jobs
  foreach my $epigenome_id (@all_epigenome_ids) {
  
    my @all_feature_type_ids = keys %{$unique{$epigenome_id}};
    foreach my $feature_type_id (@all_feature_type_ids) {
    
      #print "$epigenome_id - $feature_type_id\n";
      
      my $analysis_id = $unique{$epigenome_id}{$feature_type_id};
      
      my $epigenome    = $epigenome_adaptor->fetch_by_dbID($epigenome_id);
      my $feature_type = $feature_type_adaptor->fetch_by_dbID($feature_type_id);
      my $analysis     = $analysis_adaptor->fetch_by_dbID($analysis_id);
      
      $self->dataflow_output_id(
        {
          'species'                   => $species,
          'epigenome_id'              => $epigenome_id,
          'feature_type_id'           => $feature_type_id,
          'epigenome_production_name' => $epigenome->production_name,
          'feature_type_name'         => $feature_type->name,
          'analysis_logic_name'       => $analysis->logic_name,
          'assembly'                  => $assembly,
          'foo'                  => 'bar',
        }, 
        BRANCH_JOB_OUTPUT
      );
    }
  }
}

1;
