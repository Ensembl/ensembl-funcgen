package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::AlignFastqFile;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;

use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw( :all );

sub run {

  my $self = shift;
  
  my $tempdir          = $self->param_required('tempdir');
  my $bam_file         = $self->param_required('bam_file');
  my $species          = $self->param_required('species');
  my $to_gender        = $self->param_required('to_gender');
  my $to_assembly      = $self->param_required('to_assembly');
  my $ensembl_analysis = $self->param_required('ensembl_analysis');
  my $fastq_file       = $self->param_required('fastq_file');
  
  use Bio::EnsEMBL::Registry;
  Bio::EnsEMBL::Registry->set_disconnect_when_inactive;
  
  use Bio::EnsEMBL::Funcgen::Utils::RefBuildFileLocator;
  my $bwa_index_locator = Bio::EnsEMBL::Funcgen::Utils::RefBuildFileLocator->new;
  
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
  
  my $align_commands;
  my $number_of_fastq_files = @$fastq_file;
  
  if ($ensembl_analysis eq ENSEMBL_HODGEPODGE_ALIGNMENT_ANALYSIS) {
  
    # This means it could be either single- or paired end. The decision is 
    # made here on a case by case basis. It really shouldn't be happening at 
    # all.
    #
    my $number_of_read_files = @$fastq_file;
    
    if ($number_of_read_files == 1) {
      $ensembl_analysis = ENSEMBL_SINGLE_END_ALIGNMENT_ANALYSIS
    }
    if ($number_of_read_files == 2) {
      $ensembl_analysis = ENSEMBL_PAIRED_END_ALIGNMENT_ANALYSIS
    }
  }
  if ($ensembl_analysis eq ENSEMBL_SINGLE_END_ALIGNMENT_ANALYSIS) {
  
    my $expected_number_of_fastq_files = 1;
    if ($number_of_fastq_files != $expected_number_of_fastq_files) {
      $self->throw("Wrong number of fastq files! Got $number_of_fastq_files, but expected $expected_number_of_fastq_files.");
    }
    $align_commands = $self->create_samse_commands({
      tempdir              => $tempdir,
      samtools_fasta_index => $samtools_fasta_index,
      sorted_file          => $bam_file,
      fastq_file           => $fastq_file->[0],
      bwa_index            => $bwa_index,
    });
    
  }
  if ($ensembl_analysis eq ENSEMBL_PAIRED_END_ALIGNMENT_ANALYSIS) {
     
    my $expected_number_of_fastq_files = 2;
    if ($number_of_fastq_files != $expected_number_of_fastq_files) {
      $self->throw("Wrong number of fastq files! Got $number_of_fastq_files, but expected $expected_number_of_fastq_files.");
    }
       
    $align_commands = $self->create_sampe_commands({
      tempdir              => $tempdir,
      samtools_fasta_index => $samtools_fasta_index,
      sorted_file          => $bam_file,
      fastq_1_file         => $fastq_file->[0],
      fastq_2_file         => $fastq_file->[1],
      bwa_index            => $bwa_index,
    });
  }
  
  if (! $align_commands) {
    $self->throw(
      "Don't know how to create commands for ensembl analysis " 
      . $ensembl_analysis
    );
  }
  
  my $run_options = {
    use_bash_pipefail => 1
  };

  $self->say_with_header(Dumper($align_commands));
  
  foreach my $cmd (@$align_commands) {
  
    my $has_failed = $self->run_system_command($cmd, $run_options);
    if ($has_failed) {
      $self->throw("The following command failed:\n" . $cmd)
    }
  }
  
  if (! -e $bam_file) {
    $self->throw("The bam file $bam_file should have been created, but it doesn't exist!");
  }
  
  my $bam_file_content = `samtools view $bam_file | head`;
  
  if ($bam_file_content eq "") {
    $self->throw("The bam file $bam_file is empty!");
  }

  # Give file system time to sync
  sleep(20);
  return;
}

sub create_samse_commands {

  my $self  = shift;
  my $param = shift;
  
  my $tempdir              = $param->{tempdir};
  my $samtools_fasta_index = $param->{samtools_fasta_index};
  my $sorted_file          = $param->{sorted_file};
  my $fastq_file           = $param->{fastq_file};
  my $bwa_index            = $param->{bwa_index};
  
  my $suffix_array_index = $tempdir . '/'. 'suffix_array_index.sai';
  my $samse_file         = $tempdir . '/'. 'samse.sam';
  my $unsorted_file      = $tempdir . '/'. 'unsorted.bam';
  
  my @align_commands = (
    "bwa aln   $bwa_index $fastq_file                       > $suffix_array_index",
    "bwa samse $bwa_index $suffix_array_index $fastq_file   > $samse_file",
    "samtools view -t $samtools_fasta_index -bh $samse_file > $unsorted_file",
    "samtools sort $unsorted_file -o $sorted_file",
    "check_bam_file_has_EOF_marker.pl --bam_file $sorted_file"
  );
  return \@align_commands;
}

sub create_sampe_commands {

  my $self  = shift;
  my $param = shift;
  
  my $tempdir              = $param->{tempdir};
  my $samtools_fasta_index = $param->{samtools_fasta_index};
  my $sorted_file          = $param->{sorted_file};
  my $fastq_1_file         = $param->{fastq_1_file};
  my $fastq_2_file         = $param->{fastq_2_file};
  my $bwa_index            = $param->{bwa_index};

  my $suffix_array_index_1 = $tempdir . '/'. 'suffix_array_index_1.sai';
  my $suffix_array_index_2 = $tempdir . '/'. 'suffix_array_index_2.sai';
  my $sampe_file           = $tempdir . '/'. 'aln-pe.sam';
  my $unsorted_file        = $tempdir . '/'. 'unsorted.bam';
  
  my @align_commands = (
    "bwa aln   $bwa_index $fastq_1_file                     > $suffix_array_index_1",
    "bwa aln   $bwa_index $fastq_2_file                     > $suffix_array_index_2",
    "bwa sampe $bwa_index $suffix_array_index_1 $suffix_array_index_2 $fastq_1_file $fastq_2_file  > $sampe_file",
    "samtools view -t $samtools_fasta_index -bh $sampe_file > $unsorted_file",
    "samtools sort $unsorted_file -o $sorted_file",
    "check_bam_file_has_EOF_marker.pl --bam_file $sorted_file"
  );
  return \@align_commands;
}

1;
