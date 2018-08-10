package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::SplitExperimentReads;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::ExecutionPlanUtils qw (
    lock_execution_plan
);

use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw ( :all );

use constant {
  BRANCH_MERGE       => 2,
  BRANCH_SPLIT_FASTQ => 3,
};

sub run {

  my $self = shift;
  my $species       = $self->param_required('species');
  my $plan          = $self->param_required('execution_plan');
  my $tempdir       = $self->param_required('tempdir');
  my $data_root_dir = $self->param_required('data_root_dir');
  
  print Dumper($plan);
  
  lock_execution_plan($plan);
  
  my $bam_file_no_duplicates = $plan
    ->{output}
    ->{real}
  ;
  my $full_path_to_deduplicated_bam = $data_root_dir . '/' . $bam_file_no_duplicates;

  my $align_plan = $plan->{input};

  my $bam_file       = $align_plan->{output}->{real};
  my $read_files     = $align_plan->{input}->{read_files};
  my $alignment_name = $align_plan->{name};
  my $bam_file_real  = $align_plan->{output}->{real};
  
  my $full_path_to_merged_bam = $data_root_dir . '/' . $bam_file_real;

  my @chunks_to_be_merged;
  my $number = 0;
  foreach my $read_file (@$read_files) {
  
    my $chunk_bam_file = $tempdir  . '/' . $alignment_name . '/fastq_bams/' .  $bam_file . '.' . $number . '.bam';
    push @chunks_to_be_merged, $chunk_bam_file;
    $number++;
  
    $self->dataflow_output_id(
        {
            'species'        => $species,
            'tempdir'        => $tempdir,
            'execution_plan' => $plan,
            'read_file'      => $read_file,
            'merged_bam'     => $chunk_bam_file,
        },
        BRANCH_SPLIT_FASTQ
    );
  };
  
  $self->dataflow_output_id(
    {
      'species'          => $species,
      'chunks'           => \@chunks_to_be_merged,
      'merged_bam'       => $full_path_to_merged_bam,
      'execution_plan'   => $plan,
      'deduplicated_bam' => $full_path_to_deduplicated_bam,
    }, 
    BRANCH_MERGE
  );
  
  return;
}

1;
