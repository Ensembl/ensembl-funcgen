package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::CreateExecutionPlan;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;

sub run {

  my $self = shift;
  
  my $species                 = $self->param_required('species');
  my $experiment_id           = $self->param_required('experiment_id');
  my $ensembl_release_version = $self->param_required('ensembl_release_version');
  my $tempdir                 = $self->param_required('tempdir');
  
  
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

  use Bio::EnsEMBL::Funcgen::PeakCallingPlan::DirectoryNameBuilder;
  my $directory_name_builder 
    = Bio::EnsEMBL::Funcgen::PeakCallingPlan::DirectoryNameBuilder
      ->new(
        -root_dir                => '',
        -species                 => $species,
        -assembly                => $default_assembly,
        -ensembl_release_version => $ensembl_release_version,
      );

  use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Director;
  my $chip_seq_analysis_director 
    = Bio::EnsEMBL::Funcgen::PeakCallingPlan::Director->new;
  
  my $execution_plan 
    = $chip_seq_analysis_director->construct_execution_plan(
      {
        species                => $species, 
        assembly               => $default_assembly, 
        experiment             => $experiment,
        directory_name_builder => $directory_name_builder
      }
    );

  use Bio::EnsEMBL::Funcgen::PeakCallingPlan::ExecutionPlanUtils qw (
        lock_execution_plan
        resolve_nonterminal_symbols
        summarise
  );

  my $execution_plan_expanded = resolve_nonterminal_symbols($execution_plan);
  
#   use YAML qw( Dump );
# 
#   local $YAML::Indent     = 8;
#   local $YAML::UseAliases = 0;

  $Data::Dumper::Deepcopy = 1;
  $Data::Dumper::Sortkeys = 1;

  use Bio::EnsEMBL::Funcgen::ExecutionPlan;
  my $execution_plan_obj = Bio::EnsEMBL::Funcgen::ExecutionPlan->new(
    -experiment_id  => $experiment_id,
    -execution_plan => Dumper($execution_plan_expanded),
    -time           => time
  );

  my $execution_plan_adaptor = Bio::EnsEMBL::Registry->get_adaptor(
    $species, 
    'funcgen', 
    'ExecutionPlan'
  );
  $execution_plan_adaptor->store($execution_plan_obj);
  
  # Execution plans are written to disk to allow easier debugging.
  #
  my $dump_directory = join '/', 
    $tempdir,
    'execution_plans',
    $species,
    $experiment->get_Epigenome->production_name,
  ;
  my $file_basename = $experiment->name . '.pl';
  
  my $full_file_name = $dump_directory . '/' . $file_basename;
    
  use File::Path qw( make_path );
  make_path($dump_directory);
  
  $Data::Dumper::Deepcopy = 1;
  $Data::Dumper::Sortkeys = 1;
  
  open my $fh, '>', $full_file_name || confess("Couldn't open $full_file_name!");
  $fh->print(Dumper($execution_plan_expanded));
  $fh->close;
  
  $self->say_with_header("Execution plan written to $full_file_name", 1);

  $self->dataflow_output_id( {
    'plan'    => $execution_plan,
    'species' => $species,
  }, 2);
}

1;
