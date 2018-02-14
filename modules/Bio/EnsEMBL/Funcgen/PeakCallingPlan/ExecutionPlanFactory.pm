package Bio::EnsEMBL::Funcgen::PeakCallingPlan::ExecutionPlanFactory;

use strict;
use Data::Dumper;

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_get_or_set
);

sub _constructor_parameters {
  return {
    root_dir                => 'root_dir',
    species                 => 'species',
    ensembl_release_version => 'ensembl_release_version',
  };
}

use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_get_or_set
);

sub root_dir                   { return shift->_generic_get_or_set('root_dir',                   @_); }
sub species                    { return shift->_generic_get_or_set('species',                    @_); }
sub ensembl_release_version    { return shift->_generic_get_or_set('ensembl_release_version',    @_); }
sub chip_seq_analysis_director { return shift->_generic_get_or_set('chip_seq_analysis_director', @_); }
sub default_assembly           { return shift->_generic_get_or_set('default_assembly',           @_); }
sub directory_name_builder     { return shift->_generic_get_or_set('directory_name_builder',     @_); }

sub init {

  my $self = shift;
  
  my $species = $self->species;
  
  my $experiment_adaptor  = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'Experiment');
  my $coordsystem_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'coordsystem');

  my $default_chromosome_coordsystem = $coordsystem_adaptor->fetch_by_name('chromosome');
  my $default_assembly = $default_chromosome_coordsystem->version;

  use Bio::EnsEMBL::Funcgen::PeakCallingPlan::DirectoryNameBuilder;
  my $directory_name_builder 
    = Bio::EnsEMBL::Funcgen::PeakCallingPlan::DirectoryNameBuilder
      ->new(
        -root_dir                => $self->root_dir,
        -species                 => $self->species,
        -assembly                => $default_assembly,
        -ensembl_release_version => $self->ensembl_release_version,
      );

  use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Director;
  my $chip_seq_analysis_director = Bio::EnsEMBL::Funcgen::PeakCallingPlan::Director->new;
  
  $self->chip_seq_analysis_director($chip_seq_analysis_director);
  $self->default_assembly($default_assembly);
  $self->directory_name_builder($directory_name_builder);
  
  return;
}

sub create_execution_plan_for_experiment {

  my $self       = shift;
  my $experiment = shift;
  
  my $chip_seq_analysis_director = $self->chip_seq_analysis_director;
  
  my $execution_plan 
    = $chip_seq_analysis_director->construct_execution_plan(
      {
        species                => $self->species, 
        assembly               => $self->default_assembly, 
        experiment             => $experiment,
        directory_name_builder => $self->directory_name_builder
      }
    );
  return $execution_plan;
}

1;
