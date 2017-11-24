package Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::AlignAllPlanBuilder;

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::Constants qw ( :all );

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

  my $read_file_experimental_configurations
    = $read_file_experimental_configuration_adaptor
      ->fetch_all_by_Experiment($experiment);
  
  my @names_of_reads_to_merge;
  foreach my $read_file_experimental_configuration (@$read_file_experimental_configurations) {
  
    my $read_file = $read_file_experimental_configuration->get_ReadFile;
    push @names_of_reads_to_merge, $read_file->name,
  
  }
  
  use Bio::EnsEMBL::Funcgen::Hive::RefBuildFileLocator;
  my $bwa_index_locator = Bio::EnsEMBL::Funcgen::Hive::RefBuildFileLocator->new;

  my $to_gender = 
    $bwa_index_locator
        ->epigenome_gender_to_directory_species_gender(
            $experiment->epigenome->gender
        );
  
  use Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::AlignmentPlanFactory;
  my $alignment_plan_factory = Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::AlignmentPlanFactory
  ->new(
    -names_of_reads_to_merge => \@names_of_reads_to_merge,
    -description             => 'All reads from the ' . $self->experiment_type($experiment) . ' experiment',
    -name                    => $alignment_namer->base_name_with_duplicates,
    -to_gender               => $to_gender,
    -to_assembly             => $assembly,
    -analysis                => 'align',
    -from_experiment         => $experiment->name,
    -output_real             => $alignment_namer->bam_file_with_duplicates,
    -output_stored           => $alignment_namer->bam_file_with_duplicates_stored,
    -output_format           => BAM_FORMAT,
    -is_control              => $experiment->is_control,
    -has_all_reads           => 1,
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
    $experiment_type = 'control';
  } else {
    $experiment_type = 'signal';
  }
  return $experiment_type;
}

1;
