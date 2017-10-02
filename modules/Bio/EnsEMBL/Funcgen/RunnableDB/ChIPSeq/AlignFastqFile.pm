package Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::AlignFastqFile;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;

use constant {
  BRANCH_ACCUMULATOR => 2,
};

sub run {

  my $self = shift;
  
  my $fastq_file  = $self->param_required('fastq_file');
  my $tempdir     = $self->param_required('tempdir');
  my $bam_file    = $self->param_required('bam_file');
  my $species     = $self->param_required('species');
  my $to_gender   = $self->param_required('to_gender');
  my $to_assembly = $self->param_required('to_assembly');
  
  use Bio::EnsEMBL::Funcgen::Hive::RefBuildFileLocator;
  my $bwa_index_locator = Bio::EnsEMBL::Funcgen::Hive::RefBuildFileLocator->new;
  
  my $reference_file_root = $self->param('reference_data_root_dir');

  my $bwa_index_relative = $bwa_index_locator->locate({
    species          => $species,
    assembly         => $to_assembly,
    epigenome_gender => $to_gender,
    file_type        => 'bwa_index',
  });

  my $bwa_index = $reference_file_root . '/' . $bwa_index_relative;
  
  my $samtools_fasta_index_relative = $bwa_index_locator->locate({
    species          => $species,
    assembly         => $to_assembly,
    epigenome_gender => $to_gender,
    file_type        => 'samtools_fasta_index',
  });
  
  my $samtools_fasta_index = $reference_file_root . '/' . $samtools_fasta_index_relative;

  $self->say_with_header("samtools_fasta_index = $samtools_fasta_index");
  $self->say_with_header("bwa_index            = $bwa_index");
  
  my $suffix_array_index = $tempdir . '/suffix_array_index.sai';
  my $samse_file         = $tempdir . '/samse.sam';
  my $unsorted_file      = $tempdir . '/unsorted.bam';
  my $sorted_file        = $bam_file;
  
  my @align_commands = (
    "bwa aln   $bwa_index $fastq_file                       > $suffix_array_index",
    "bwa samse $bwa_index $suffix_array_index $fastq_file   > $samse_file",
    "samtools view -t $samtools_fasta_index -bh $samse_file > $unsorted_file",
    "samtools sort $unsorted_file -o $sorted_file",
  );
  
  my $run_options = {
    use_bash_pipefail => 1
  };
  
  foreach my $cmd (@align_commands) {
  
    my $has_failed = $self->run_system_command($cmd, $run_options);
    if ($has_failed) {
      $self->throw("The following command failed:\n" . $cmd)
    }
  }
  $self->dataflow_output_id(
    {
      'chunk' => $sorted_file,
    }, 
    BRANCH_ACCUMULATOR
  );
  return;
}

1;
