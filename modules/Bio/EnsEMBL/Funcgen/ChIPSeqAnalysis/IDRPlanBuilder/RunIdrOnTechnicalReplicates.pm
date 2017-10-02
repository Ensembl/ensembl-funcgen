package Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::IDRPlanBuilder::RunIdrOnTechnicalReplicates;

use strict;
use Data::Dumper;

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

sub construct {
  my $self  = shift;
  my $param = shift;
  
  my $species                = $param->{species};
  my $experiment             = $param->{experiment};
  my $directory_name_builder = $param->{directory_name_builder};
  my $assembly               = $param->{assembly};

  my $read_file_experimental_configuration_adaptor 
  = Bio::EnsEMBL::Registry->get_adaptor(
    $species, 
    'funcgen', 
    'ReadFileExperimentalConfiguration'
  );

  my $read_file_experimental_configuration_list = 
    $read_file_experimental_configuration_adaptor
      ->fetch_all_by_Experiment($experiment);

  my @align_replicates_plan;
  foreach my $read_file_experimental_configuration 
    (@$read_file_experimental_configuration_list) {
  
    my $read_file = $read_file_experimental_configuration->get_ReadFile;
    my $technical_replicate_number
      = $read_file_experimental_configuration
        ->technical_replicate;
    
    my $to_gender = $experiment->epigenome->gender;

    use Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::AlignmentNamer;
    my $alignment_namer = Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::AlignmentNamer
      ->new(
          -directory_name_builder => $directory_name_builder,
          -experiment                  => $experiment,
          -biological_replicate_number => undef,
          -technical_replicate_number  => $technical_replicate_number,
      );
    push @align_replicates_plan, {
      remove_duplicates => {
        align => {
        
          # Using array with one item to be consistent. align
          # should always point to an array of read files.
          #
          read_files  => [ $read_file->name, ],
          type        => 'individual read files',
          name        => $alignment_namer->base_name_with_duplicates,
          bam_file    => {
            real   => $alignment_namer->bam_file_with_duplicates,
            stored => $alignment_namer->bam_file_with_duplicates_stored,
          },
          to_gender   => $to_gender,
          to_assembly => $assembly,
        },
        name        => $alignment_namer->base_name_no_duplicates,
        bam_file    => {
          real   => $alignment_namer->bam_file_no_duplicates,
          stored => $alignment_namer->bam_file_no_duplicates_stored,
        },
        bigwig_file => {
          real   => $alignment_namer->bigwig_file_no_duplicates,
          stored => $alignment_namer->bigwig_file_no_duplicates_stored,
        },
      },
    };
  }
  
  my $idr_plan = {
    alignment_replicates => \@align_replicates_plan,
    type => Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::IDRStrategy
      ->RUN_IDR_ON_TECHNICAL_REPLICATES
  };
  return $idr_plan;
}

1;
