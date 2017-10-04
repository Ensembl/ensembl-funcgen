package Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::IDRPlanBuilder::RunIdrOnBiologicalReplicates;

use strict;
use Data::Dumper;

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

  my $biological_replicate_numbers = 
    $read_file_experimental_configuration_adaptor
      ->fetch_all_biological_replicate_numbers_from_Experiment($experiment);

  my @align_technical_replicates_plan;

  my $alignment = [];
  
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

    my $to_gender = $experiment->epigenome->gender;

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
      -output_real             => $alignment_namer->bam_file_with_duplicates,
      -output_stored           => $alignment_namer->bam_file_with_duplicates_stored,
      -output_format           => 'bam',
    );

    my $alignment_plan = $alignment_plan_factory->product;
    
    my $remove_duplicates_plan_factory
      = Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::RemoveDuplicatesPlanFactory
        ->new(
          -input         => create_ref($alignment_plan),
          -name          => $alignment_namer->base_name_no_duplicates,
          -output_real   => $alignment_namer->bam_file_no_duplicates,
          -output_stored => $alignment_namer->bam_file_no_duplicates_stored,
          -output_format => 'bam',
        );
    my $remove_duplicates_plan = $remove_duplicates_plan_factory->product;
    
    push @align_technical_replicates_plan, create_ref($remove_duplicates_plan);
    
    push @$alignment, $alignment_plan;
    push @$alignment, $remove_duplicates_plan;
  }
  
  my $idr_plan = {
    alignment_replicates => \@align_technical_replicates_plan,
    type => 'idr',
    name => $experiment->name,
    strategy => Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::IDRStrategy
      ->RUN_IDR_ON_BIOLOGICAL_REPLICATES
  };
  
  $self->set_Alignment($alignment);
  
  return $idr_plan;
}

sub create_ref {
  my $plan = shift;
  
  my $name = $plan->{name};
  my $type = $plan->{type};
  
  if (! defined $name) {
    die;
  }
  if (! defined $type) {
    die;
  }
  my $ref = '#' . $type . ':' . $name . '#';
  return $ref;
}

1;
