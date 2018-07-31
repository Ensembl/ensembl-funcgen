package Bio::EnsEMBL::Funcgen::PeakCallingPlan::Director;

use strict;
use Data::Dumper;

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

use Bio::EnsEMBL::Funcgen::PeakCallingPlan::IDRPlanBuilderFactory;
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::PlanMetadataBuilder;
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::AlignAllPlanBuilder;
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::PeakCallingPlanBuilder;
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::SignalFilePlanBuilder;
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::PlanMetadataBuilder;
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::RemoveDuplicatesPlanBuilder;
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw ( :all );
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::ExecutionPlanUtils qw ( create_ref );

use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_get_or_set
);

sub construct_execution_plan {

  my $self  = shift;
  my $param = shift;
  
  my $species                = $param->{species};
  my $experiment             = $param->{experiment};
  my $directory_name_builder = $param->{directory_name_builder};
  my $assembly               = $param->{assembly};

  use Bio::EnsEMBL::Funcgen::PeakCallingPlan::AlignmentNamer;
  my $alignment_namer = Bio::EnsEMBL::Funcgen::PeakCallingPlan::AlignmentNamer->new(
      -directory_name_builder      => $directory_name_builder,
      -experiment                  => $experiment,
  );
  $param->{alignment_namer} = $alignment_namer;

  my $control_experiment = $param->{experiment}->get_control;
  
  my $control_alignment_namer;
  
  if (defined $control_experiment) {
    $control_alignment_namer = Bio::EnsEMBL::Funcgen::PeakCallingPlan::AlignmentNamer->new(
        -directory_name_builder      => $directory_name_builder,
        -experiment                  => $experiment->get_control,
    );
  }

  use Bio::EnsEMBL::Funcgen::Utils::RefBuildFileLocator;
  my $bwa_index_locator = Bio::EnsEMBL::Funcgen::Utils::RefBuildFileLocator->new;
  
  my $samtools_fasta_index = $bwa_index_locator->locate({
    species          => $species,
    epigenome_gender => $experiment->epigenome->gender,
    assembly         => $assembly,
    file_type        => 'samtools_fasta_index',
  });

  my $chromosome_lengths_by_species_assembly = $bwa_index_locator->locate({
    species          => $species,
    epigenome_gender => $experiment->epigenome->gender,
    assembly         => $assembly,
    file_type        => 'chromosome_lengths_by_species_assembly',
  });

  #
  # Create metadata
  #
  
  my $plan_meta_data_builder 
    = Bio::EnsEMBL::Funcgen::PeakCallingPlan::PlanMetadataBuilder
      ->new;
  
  $plan_meta_data_builder->construct($param);
  my $meta_data = $plan_meta_data_builder->get_meta_data;
    
  #
  # Align all reads for peak calling
  #
  
  my $align_all_plan_builder 
    = Bio::EnsEMBL::Funcgen::PeakCallingPlan::AlignAllPlanBuilder
        ->new;
  
  $align_all_plan_builder->construct($param);
      
  my $align_all_read_files_for_experiment_plan 
    = $align_all_plan_builder
      ->get_Alignment;

  my %control_param = %$param;
  
  my $align_all_read_files_for_control_plan;
  
  if (defined $control_experiment) {
  
    $control_param{experiment}      = $param->{experiment}->get_control;
    $control_param{alignment_namer} = $control_alignment_namer;
    
    $align_all_plan_builder
        ->construct(\%control_param);
        
    $align_all_read_files_for_control_plan
      = $align_all_plan_builder
        ->get_Alignment;
  } else {
    $align_all_read_files_for_control_plan = {
        name        => NO_CONTROL_FLAG,
        type        => ALIGNMENT_TYPE,
        description => 'No control alignment for this experiment.',
    };
  }
  
  #
  # Remove duplicates
  #

  my $remove_duplicates_plan_builder
    = Bio::EnsEMBL::Funcgen::PeakCallingPlan::RemoveDuplicatesPlanBuilder
      ->new;
  
  $remove_duplicates_plan_builder->set_input         (create_ref($align_all_read_files_for_experiment_plan));
  $remove_duplicates_plan_builder->set_name          ($alignment_namer->base_name_no_duplicates);
  $remove_duplicates_plan_builder->set_is_control    ( FALSE );
  $remove_duplicates_plan_builder->set_output_real   ($alignment_namer->bam_file_no_duplicates);
  $remove_duplicates_plan_builder->set_output_stored ($alignment_namer->bam_file_no_duplicates_stored);
  $remove_duplicates_plan_builder->set_output_format ( BAM_FORMAT );
  $remove_duplicates_plan_builder->set_experiment    ( $experiment->name );
  
  $remove_duplicates_plan_builder->construct;
  
  my $remove_duplicates_plan = $remove_duplicates_plan_builder->get_plan;
  
  my $remove_duplicates_from_control_plan;
  if (defined $control_experiment) {
  
    $remove_duplicates_plan_builder->set_input         (create_ref($align_all_read_files_for_control_plan));
    $remove_duplicates_plan_builder->set_name          ($control_alignment_namer->base_name_no_duplicates);
    $remove_duplicates_plan_builder->set_output_real   ($control_alignment_namer->bam_file_no_duplicates);
    $remove_duplicates_plan_builder->set_output_stored ($control_alignment_namer->bam_file_no_duplicates_stored);
    $remove_duplicates_plan_builder->set_is_control    ( TRUE );
    $remove_duplicates_plan_builder->set_output_format ( BAM_FORMAT );
    $remove_duplicates_plan_builder->set_experiment    ( $experiment->get_control->name );

    $remove_duplicates_plan_builder->construct;
    
    $remove_duplicates_from_control_plan = $remove_duplicates_plan_builder->get_plan;
  } else {
    $remove_duplicates_from_control_plan = {
        name        => NO_CONTROL_FLAG,
        type        => REMOVE_DUPLICATES_ANALYSIS,
        description => 'No control alignment for this experiment.',
    };
  }
  
  #
  # Bigwig files
  #
  
  my $signal_file_plan_builder
    = Bio::EnsEMBL::Funcgen::PeakCallingPlan::SignalFilePlanBuilder
      ->new;
  
  $signal_file_plan_builder->set_alignment        (create_ref($remove_duplicates_plan));
  $signal_file_plan_builder->set_alignment_namer  ($alignment_namer);
  $signal_file_plan_builder->set_is_control ( FALSE );
  
  $signal_file_plan_builder->construct;
  
  my $signal_file_plan = $signal_file_plan_builder->get_signal_plan;
  
  my $control_file_plan;
  if (defined $control_experiment) {
    $signal_file_plan_builder->set_alignment        (create_ref($remove_duplicates_from_control_plan));
    $signal_file_plan_builder->set_alignment_namer  ($control_alignment_namer);
    $signal_file_plan_builder->set_is_control ( TRUE );
    
    $signal_file_plan_builder->construct;
    
    $control_file_plan = $signal_file_plan_builder->get_signal_plan;
  } else {
    $control_file_plan = {
        name        => NO_CONTROL_FLAG,
        type        => SIGNAL_EXPERIMENT,
        description => 'No control alignment for this experiment.',
    };
  }
  
  #
  # IDR plan
  #
  
  my $idr_plan_builder_factory
    = Bio::EnsEMBL::Funcgen::PeakCallingPlan::IDRPlanBuilderFactory
      ->new;
  
  my $idr_plan_builder 
    = $idr_plan_builder_factory
      ->make($experiment);

  $idr_plan_builder->construct($param);
  
  my $idr_plan           = $idr_plan_builder->get_idr_plan;
  my $alignments_for_idr = $idr_plan_builder->get_Alignments;

  #
  # Peak calling plan
  #

  my $peak_calling_plan_builder 
    = Bio::EnsEMBL::Funcgen::PeakCallingPlan::PeakCallingPlanBuilder
      ->new;
  
  my $peaks_output_dir 
    = $directory_name_builder
      ->peak_calling_output_dir_by_Experiment(
        $experiment
      );
  
  $peak_calling_plan_builder->set_signal_alignment        ( create_ref($remove_duplicates_plan) );
  $peak_calling_plan_builder->set_control_alignment       ( create_ref($remove_duplicates_from_control_plan) );
  $peak_calling_plan_builder->set_idr                     ( $idr_plan );
  $peak_calling_plan_builder->set_peaks_output_dir        ( $peaks_output_dir );
  $peak_calling_plan_builder->set_alignment_namer         ( $alignment_namer );
  if ($control_alignment_namer) {
    $peak_calling_plan_builder->set_control_alignment_namer ( $control_alignment_namer );
  }
  $peak_calling_plan_builder->set_experiment              ( $experiment );
  $peak_calling_plan_builder->set_samtools_fasta_index    ( $samtools_fasta_index );
  $peak_calling_plan_builder
    ->set_chromosome_lengths_by_species_assembly
      ( $chromosome_lengths_by_species_assembly );
  
  
  $peak_calling_plan_builder->construct;
  
  my $peak_calling_plan    = $peak_calling_plan_builder->get_peakcalling_plan;
  my $bam_file_to_bed_plan = $peak_calling_plan_builder->get_bam_to_bed_conversion_plan;
  
  # Assemble components into final execution plan
  
  my @all_alignment_plans = (
    $align_all_read_files_for_experiment_plan,
    $remove_duplicates_plan,
    $align_all_read_files_for_control_plan,
    $remove_duplicates_from_control_plan,
    @$alignments_for_idr
  );
  
  my $alignment_plan = {};
  
  foreach my $current_alignment_plan (@all_alignment_plans) {
    $alignment_plan->{$current_alignment_plan->{name}} = $current_alignment_plan;
  }
  
  my $bigwig_plan = {
    $signal_file_plan ->{name} => $signal_file_plan,
    $control_file_plan->{name} => $control_file_plan,
  };
  
  my $execution_plan = {
    meta_data  => $meta_data,
    signal     => $bigwig_plan,
    call_peaks => $peak_calling_plan,
    alignment  => $alignment_plan,
    idr        => $idr_plan,
    bam_to_bed => $bam_file_to_bed_plan
  };
  return $execution_plan;
}

1;
