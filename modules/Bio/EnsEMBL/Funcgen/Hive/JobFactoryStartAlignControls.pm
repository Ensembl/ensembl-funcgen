package Bio::EnsEMBL::Funcgen::Hive::JobFactoryStartAlignControls;

use warnings;
use strict;

sub fetch_input {  
  my $self = shift;
  
  $self->SUPER::fetch_input();
  my $result_set = $self->fetch_Set_input('ResultSet');

  my $run_controls = $self->get_param_method('result_set_groups', 'silent') ? 1 : 0;
  $self->set_param_method('run_controls', $run_controls);
  my $merge        = $self->get_param_method('merge', 'silent', $run_controls); 
  $self->get_param_method('checksum_optional', 'silent');

  $self->get_param_method('fastq_chunk_size', 'silent', '16000000');#default should run in 30min-1h 
  
  $self->get_output_work_dir_methods($self->alignment_dir($result_set, 1, $run_controls));#default output_dir 
  return;
}

sub run {
  my $self         = shift;
  my $result_set   = $self->ResultSet;
  my $run_controls = $self->run_controls;
  my $merge        = $self->merge;
 
  my %batch_params = %{$self->batch_params};
  
  # Create funnel job first so hive can start creating the jobs in the database.
  #
  my %signal_info;  
  if($run_controls) {
    # for flow to MergeControlAlignments_and_QC
    %signal_info = (result_set_groups => $self->result_set_groups);
  }
  
  my $file_prefix  = $self->get_alignment_path_prefix_by_ResultSet($result_set, $run_controls); 
  
  my $bam_file_with_unmapped_reads_and_duplicates = $file_prefix.'.with_unmapped_reads_and_duplicates.bam';
  
  # Name of the final deduplicated bam file
  #
  my $bam_file = $file_prefix . '.bam';

  $self->dataflow_output_id({
    %batch_params,
    set_type   => 'ResultSet',
    set_name   => $result_set->name,
    dbID       => $result_set->dbID,
    fastq_files => $self->fastq_files,
    output_dir => $self->output_dir,
    set_prefix => $set_prefix,
    bam_file_with_unmapped_reads_and_duplicates => $bam_file_with_unmapped_reads_and_duplicates,
    bam_file => $bam_file,
	%signal_info
      }, 3
  );

  return;
}

1;
