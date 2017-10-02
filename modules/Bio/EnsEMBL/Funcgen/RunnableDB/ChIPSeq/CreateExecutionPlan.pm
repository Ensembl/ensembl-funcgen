package Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::CreateExecutionPlan;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;

sub run {

  my $self = shift;
  
  my $species                 = $self->param_required('species');
  my $experiment_id           = $self->param_required('experiment_id');
  my $ensembl_release_version = $self->param_required('ensembl_release_version');
  
  my $experiment_adaptor = Bio::EnsEMBL::Registry->get_adaptor(
    $species, 
    'funcgen', 
    'Experiment'
  );
  my $coordsystem_adaptor = Bio::EnsEMBL::Registry->get_adaptor(
    $species, 
    'core', 
    'coordsystem'
  );
  my $experiment = $experiment_adaptor
    ->fetch_by_dbID($experiment_id);

  my $default_chromosome_coordsystem = $coordsystem_adaptor
    ->fetch_by_name('chromosome');
    
  my $default_assembly = $default_chromosome_coordsystem
    ->version;

  use Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::DirectoryNameBuilder;
  my $directory_name_builder 
    = Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::DirectoryNameBuilder
      ->new(
        -root_dir                => '',
        -species                 => $species,
        -assembly                => $default_assembly,
        -ensembl_release_version => $ensembl_release_version,
      );

  use Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::Director;
  my $chip_seq_analysis_director 
    = Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::Director->new;
  
  my $execution_plan 
    = $chip_seq_analysis_director->construct_execution_plan(
      {
        species                => $species, 
        assembly               => $default_assembly, 
        experiment             => $experiment,
        directory_name_builder => $directory_name_builder
      }
    );
  
  $self->dataflow_output_id( {
    'plan'    => $execution_plan,
    'species' => $species,
  }, 2);
}

1;
