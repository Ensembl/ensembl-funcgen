package Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::IDRPlanBuilder::RunIdrOnBiologicalReplicates;

use strict;
use Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::ExecutionPlanUtils qw ( create_ref );
use Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::Constants qw( :all );
use Data::Dumper;

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

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
  my $self = shift;
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

  my $biological_replicate_numbers = 
    $read_file_experimental_configuration_adaptor
      ->fetch_all_biological_replicate_numbers_from_Experiment($experiment);

  my @align_technical_replicates_plan;

  my $alignment = [];

  my $remove_duplicates_plan_builder
    = Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::RemoveDuplicatesPlanBuilder
      ->new;

  $remove_duplicates_plan_builder->set_output_format ( BAM_FORMAT );
  $remove_duplicates_plan_builder->set_experiment    ( $experiment->name );

  foreach my $biological_replicate_number (@$biological_replicate_numbers) {
  
    my $read_file_experimental_configurations = 
      $read_file_experimental_configuration_adaptor
      ->fetch_all_technical_replicates_by_Experiment_and_biological_replicate_number(
        $experiment, 
        $biological_replicate_number
    );
    
    my @merge_technical_replicates;
    foreach my $read_file_experimental_configuration 
      (@$read_file_experimental_configurations) {
      
      my $read_file = $read_file_experimental_configuration->get_ReadFile;
      push @merge_technical_replicates, $read_file->name;
    }

    use Bio::EnsEMBL::Funcgen::Hive::RefBuildFileLocator;
    my $bwa_index_locator = Bio::EnsEMBL::Funcgen::Hive::RefBuildFileLocator->new;

    my $to_gender = 
        $bwa_index_locator
            ->epigenome_gender_to_directory_species_gender(
                $experiment->epigenome->gender
            );

    $alignment_namer->biological_replicate_number($biological_replicate_number);
    $alignment_namer->technical_replicate_number(undef);
    
    use Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::AlignmentPlanFactory;
    my $alignment_plan_factory = Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::AlignmentPlanFactory
    ->new(
      -names_of_reads_to_merge => \@merge_technical_replicates,
      -description             => "Biological replicate $biological_replicate_number",
      -name                    => $alignment_namer->base_name_with_duplicates,
      -to_gender               => $to_gender,
      -to_assembly             => $assembly,
      -analysis                => 'align',
      -from_experiment         => $experiment->name,
      -output_real             => $alignment_namer->bam_file_with_duplicates,
      -output_stored           => $alignment_namer->bam_file_with_duplicates_stored,
      -output_format           => BAM_FORMAT,
      -is_control              => $experiment->is_control,
      -has_all_reads           => 0,
    );

    my $alignment_plan = $alignment_plan_factory->product;
    
    $remove_duplicates_plan_builder->set_input         (create_ref($alignment_plan));
    $remove_duplicates_plan_builder->set_name          ($alignment_namer->base_name_no_duplicates,);
    $remove_duplicates_plan_builder->set_output_real   ($alignment_namer->bam_file_no_duplicates);
    $remove_duplicates_plan_builder->set_output_stored ($alignment_namer->bam_file_no_duplicates_stored);
    
    $remove_duplicates_plan_builder->construct;
    
    my $remove_duplicates_plan = $remove_duplicates_plan_builder->get_plan;
    
    push @align_technical_replicates_plan, create_ref($remove_duplicates_plan);
    
    push @$alignment, $alignment_plan;
    push @$alignment, $remove_duplicates_plan;
  }
  
  my $idr_plan = {
    alignment_replicates => \@align_technical_replicates_plan,
    type                 => 'idr',
    name                 => $experiment->name,
    strategy             => RUN_IDR_ON_BIOLOGICAL_REPLICATES
  };
  
  $self->set_Alignments($alignment);
  $self->set_idr_plan($idr_plan);
  
  return;
}

1;
