package Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::AlignAllPlanBuilder;

use strict;
use Data::Dumper;
use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

sub construct {
  my $self = shift;
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

  my $read_file_experimental_configurations
    = $read_file_experimental_configuration_adaptor
      ->fetch_all_by_Experiment($experiment);
  
  my @names_of_reads_to_merge;
  foreach my $read_file_experimental_configuration (@$read_file_experimental_configurations) {
  
    my $read_file = $read_file_experimental_configuration->get_ReadFile;
    push @names_of_reads_to_merge, $read_file->name,
  
  }
  my $to_gender = $experiment->epigenome->gender;
  
  use Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::AlignmentNamer;
  my $alignment_namer = Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::AlignmentNamer->new(
      -directory_name_builder      => $directory_name_builder,
      -experiment                  => $experiment,
      -biological_replicate_number => undef,
      -technical_replicate_number  => undef,
  );
  my $align_plan = {
      remove_duplicates => {
        align => {
          read_files => \@names_of_reads_to_merge,
          type => 'all reads from the experiment',
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
    }
  };
  return $align_plan;
}

1;
