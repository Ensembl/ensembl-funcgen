package Bio::EnsEMBL::Funcgen::PeakCallingPlan::IDRPlanBuilder::RunIdrOnTechnicalReplicates;

use strict;
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw( :all );
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::ExecutionPlanUtils qw ( create_ref );
use Data::Dumper;

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';
with 'Bio::EnsEMBL::Funcgen::PeakCallingPlan::select_EnsemblAlignmentAnalysis';
with 'Bio::EnsEMBL::Funcgen::PeakCallingPlan::summarise_ReadFile';

use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_set
  _generic_get
);

sub set_Alignments {
  my $self = shift;
  my $obj  = shift;
  
  return $self->_generic_set('alignment', undef, $obj);
}

sub get_Alignments { 
  return shift->_generic_get('alignment');
}

sub set_idr_plan {
  my $self = shift;
  my $obj  = shift;
  return $self->_generic_set('idr_plan', undef, $obj);
}

sub get_idr_plan { 
  return shift->_generic_get('idr_plan');
}

sub construct {
  my $self  = shift;
  my $param = shift;
  
  my $species                = $param->{species};
  my $experiment             = $param->{experiment};
  my $assembly               = $param->{assembly};
  my $alignment_namer        = $param->{alignment_namer};

  my $read_file_experimental_configuration_adaptor 
  = Bio::EnsEMBL::Registry->get_adaptor(
    $species, 
    'funcgen', 
    'ReadFileExperimentalConfiguration'
  );

  my $read_file_adaptor 
  = Bio::EnsEMBL::Registry->get_adaptor(
    $species, 
    'funcgen', 
    'ReadFile'
  );

  my $read_file_experimental_configuration_list = 
    $read_file_experimental_configuration_adaptor
      ->fetch_all_by_Experiment($experiment);

  my $alignment = [];

  my $remove_duplicates_plan_builder
    = Bio::EnsEMBL::Funcgen::PeakCallingPlan::RemoveDuplicatesPlanBuilder
      ->new;
  $remove_duplicates_plan_builder->set_output_format ( BAM_FORMAT );
  $remove_duplicates_plan_builder->set_experiment    ( $experiment->name );

  my $has_single_end_reads_only = 1;
  my $has_paired_end_reads_only = 1;
  
  READ_FILE_EXPERIMENTAL_CONFIGURATION:
  foreach my $read_file_experimental_configuration 
    (@$read_file_experimental_configuration_list) {
    
    my $read_file = $read_file_experimental_configuration->get_ReadFile;
    
    if ($read_file->is_paired_end) {
      $has_single_end_reads_only = undef;
    }
    if (! $read_file->is_paired_end) {
      $has_paired_end_reads_only = undef;
    }
  }

  my $ensembl_alignment_analysis = select_EnsemblAlignmentAnalysis({
    experiment => $experiment,
    has_single_end_reads_only => $has_single_end_reads_only,
    has_paired_end_reads_only => $has_paired_end_reads_only,
  });

  my @align_replicates_plan;
  foreach my $read_file_experimental_configuration 
    (@$read_file_experimental_configuration_list) {
  
    my $read_file = $read_file_experimental_configuration->get_ReadFile;
    my $technical_replicate_number
      = $read_file_experimental_configuration
        ->technical_replicate;
    
    use Bio::EnsEMBL::Funcgen::Utils::RefBuildFileLocator;
    my $bwa_index_locator = Bio::EnsEMBL::Funcgen::Utils::RefBuildFileLocator->new;

    my $to_gender = 
        $bwa_index_locator
            ->epigenome_gender_to_directory_species_gender(
                $experiment->epigenome->gender
            );

    $alignment_namer->unset_biological_replicate_number;
    $alignment_namer->technical_replicate_number($technical_replicate_number);

    use Bio::EnsEMBL::Funcgen::PeakCallingPlan::AlignmentPlanFactory;
    my $alignment_plan_factory = Bio::EnsEMBL::Funcgen::PeakCallingPlan::AlignmentPlanFactory
    ->new(
      -names_of_reads_to_merge => [ 
        summarise_ReadFile( 
          $read_file_adaptor->fetch_by_name(
            $read_file->name
          )
        )
      ],
      -description             => "Technical replicate $technical_replicate_number",
      -name                    => $alignment_namer->base_name_with_duplicates,
      -to_gender               => $to_gender,
      -to_assembly             => $assembly,
      -ensembl_analysis        => $ensembl_alignment_analysis,
      -from_experiment         => $experiment->name,
      -output_real             => $alignment_namer->bam_file_with_duplicates,
      -output_stored           => $alignment_namer->bam_file_with_duplicates_stored,
      -output_format           => BAM_FORMAT,
      -is_control              => $experiment->is_control,
      -has_all_reads           => FALSE,
    );
    
    my $alignment_plan = $alignment_plan_factory->product;
    
    $remove_duplicates_plan_builder->set_input         (create_ref($alignment_plan));
    $remove_duplicates_plan_builder->set_name          ($alignment_namer->base_name_no_duplicates,);
    $remove_duplicates_plan_builder->set_output_real   ($alignment_namer->bam_file_no_duplicates);
    $remove_duplicates_plan_builder->set_output_stored ($alignment_namer->bam_file_no_duplicates_stored);
    
    $remove_duplicates_plan_builder->construct;
    
    my $remove_duplicates_plan = $remove_duplicates_plan_builder->get_plan;
    
    push @align_replicates_plan, create_ref($remove_duplicates_plan);
    
    push @$alignment, $alignment_plan;
    push @$alignment, $remove_duplicates_plan;
  }
  
  my $idr_plan = {
    alignment_replicates => \@align_replicates_plan,
    type                 => IDR_ANALYSIS,
    name                 => $experiment->name,
    strategy             => RUN_IDR_ON_TECHNICAL_REPLICATES
  };
  
  $self->set_Alignments($alignment);
  $self->set_idr_plan($idr_plan);
  
  return $idr_plan;
}

1;

