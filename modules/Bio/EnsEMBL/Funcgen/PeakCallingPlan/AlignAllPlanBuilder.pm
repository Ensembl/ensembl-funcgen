package Bio::EnsEMBL::Funcgen::PeakCallingPlan::AlignAllPlanBuilder;

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw ( :all );

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_set
  _generic_get
);

sub set_Alignment {
  my $self = shift;
  my $obj  = shift;
  
  return $self->_generic_set('alignment', undef, $obj);
}

sub get_Alignment { 
  return shift->_generic_get('alignment');
}

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::PeakCallingPlan::summarise_ReadFile';
with 'Bio::EnsEMBL::Funcgen::PeakCallingPlan::select_EnsemblAlignmentAnalysis';

sub construct {
  my $self = shift;
  my $param = shift;
  
  my $species                = $param->{species};
  my $experiment             = $param->{experiment};
  my $directory_name_builder = $param->{directory_name_builder};
  my $assembly               = $param->{assembly};
  my $alignment_namer        = $param->{alignment_namer};
  
  my $read_file_experimental_configuration_adaptor 
  = Bio::EnsEMBL::Registry->get_adaptor(
    $species, 
    'funcgen', 
    'ReadFileExperimentalConfiguration'
  );

  my $read_file_experimental_configurations
    = $read_file_experimental_configuration_adaptor
      ->fetch_all_by_Experiment($experiment);
  
  my $has_single_end_reads_only = 1;
  my $has_paired_end_reads_only = 1;
  
  my @names_of_reads_to_merge;
  READ_FILE_EXPERIMENTAL_CONFIGURATION:
  foreach my $read_file_experimental_configuration (@$read_file_experimental_configurations) {
  
    my $read_file = $read_file_experimental_configuration->get_ReadFile;
    
    if ($read_file->is_paired_end) {
      $has_single_end_reads_only = undef;
    }
    if (! $read_file->is_paired_end) {
      $has_paired_end_reads_only = undef;
    }
    
    if ($read_file->is_paired_end && $read_file->paired_end_tag != 1) {
      next READ_FILE_EXPERIMENTAL_CONFIGURATION;
    }
    
    push @names_of_reads_to_merge, summarise_ReadFile($read_file);
    next READ_FILE_EXPERIMENTAL_CONFIGURATION;
  }
  
  my $ensembl_alignment_analysis = select_EnsemblAlignmentAnalysis({
    experiment => $experiment,
    has_single_end_reads_only => $has_single_end_reads_only,
    has_paired_end_reads_only => $has_paired_end_reads_only,
  });
  
  use Bio::EnsEMBL::Funcgen::Utils::RefBuildFileLocator;
  my $bwa_index_locator = Bio::EnsEMBL::Funcgen::Utils::RefBuildFileLocator->new;

  my $to_gender = 
    $bwa_index_locator
        ->epigenome_gender_to_directory_species_gender(
            $experiment->epigenome->gender
        );
  
    my $is_control;
    if ($experiment->is_control) {
        $is_control = TRUE;
    } else {
        $is_control = FALSE;
    }
  
  use Bio::EnsEMBL::Funcgen::PeakCallingPlan::AlignmentPlanFactory;
  my $alignment_plan_factory = Bio::EnsEMBL::Funcgen::PeakCallingPlan::AlignmentPlanFactory
  ->new(
    -names_of_reads_to_merge => \@names_of_reads_to_merge,
    -description             => 'All reads from the ' . $self->experiment_type($experiment) . ' experiment',
    -name                    => $alignment_namer->base_name_with_duplicates,
    -to_gender               => $to_gender,
    -to_assembly             => $assembly,
    -ensembl_analysis        => $ensembl_alignment_analysis,
    -from_experiment         => $experiment->name,
    -output_real             => $alignment_namer->bam_file_with_duplicates,
    -output_stored           => $alignment_namer->bam_file_with_duplicates_stored,
    -output_format           => BAM_FORMAT,
    -is_control              => $is_control,
    -has_all_reads           => TRUE,
  );
  my $align_plan = $alignment_plan_factory->product;
  
  $self->set_Alignment($align_plan);
  
  return;
}

sub experiment_type {
  my $self = shift;
  my $experiment = shift;

  my $experiment_type;
  if ($experiment->is_control) {
    $experiment_type = CONTROL_EXPERIMENT;
  } else {
    $experiment_type = SIGNAL_EXPERIMENT;
  }
  return $experiment_type;
}

1;
