
=pod 

=head1 NAME

Bio::EnsEMBL::Hive::RunnableDB::Funcgen::RunWiggleTools

=head1 DESCRIPTION
=cut

package Bio::EnsEMBL::Funcgen::Hive::RunWiggleTools;

use warnings;
use strict;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd
                                               run_backtick_cmd );
use Bio::EnsEMBL::Funcgen::Sequencing::SeqTools qw( write_chr_length_file);
use Bio::EnsEMBL::Utils::Exception         qw( throw );

use base ('Bio::EnsEMBL::Funcgen::Hive::BaseDB');

sub fetch_input {
  my $self = shift;
  $self->SUPER::fetch_input;

  if($self->get_param_method('set_type', 'silent')) {
    my $result_set = $self->fetch_Set_input('ResultSet');
    $self->helper->debug(1, "RunWiggleTools::fetch_input got ResultSet:\t".$result_set->name);

    if($self->get_param_method('input_files', 'silent')) {
      $self->throw_no_retry("It is unsafe to specify the input_files are set parameters.\n".
      'Please remove input_files from input_id');
    }
    my $is_control;
    if (
      defined $self->param('type')
      && $self->param('type') eq 'control'
    ) {
      $is_control = 1;
    }
    $self->input_files([$self->get_alignment_files_by_ResultSet_formats($result_set, ['bam'], $is_control)->{bam}]);

    my $output_prefix = $self->bigwig_output_dir . '/' . $result_set->name;
    
    if (
      defined $self->param('type')
      && $self->param('type') eq 'replicate'
    ) {
      $output_prefix = $self->bigwig_output_dir . '/' . $self->get_file_base_for_ResultSet($result_set, $is_control);
    }
    
    warn "output_prefix = $output_prefix";
#      die();
    
    $self->get_param_method(
      'output_prefix',
      'silent',
      $output_prefix
    );
    run_system_cmd("mkdir -p " . $self->bigwig_output_dir);
  }
  return;
}


####################################################
## Core function
####################################################

sub run {
  my $self        = shift;
  
  my $chromosome_length_file = write_chr_length_file($self->slice_objects);
  my $output = $self->output_prefix.'.bw';

  # bigWigToWig can't cope with colons, so replacing with underscores
  $output =~ s/:/_/g;

  $self->throw_no_retry("bigWigToWig can't cope with colons, these should not be in the experiment names! ($output)")
    if ($output =~ /:/);

  # RUN THE COMMAND
  # Presence of index files is input format specific, so not explicitly tested/generated here
  # Always use write to avoid redirects
  my $cmd = "wiggletools write - ". $self->_build_rpkm_cmd . ' | wigToBigWig -fixedSummaries stdin ' . $chromosome_length_file . ' ' . $output;
  
  $self->helper->debug(1, "Running:\n\t".$cmd);
  run_system_cmd($cmd);

  if($self->set_type){
    my $result_set = $self->ResultSet;
    $result_set->adaptor->dbfile_data_root($self->db_output_dir);
    $result_set->dbfile_path($output);
    $result_set->adaptor->store_dbfile_path($result_set, 'BIGWIG');
  }
  return;
}

sub _build_rpkm_cmd {
  my $self = shift;

  # This assumes bam input with index
  my @one_bam = @{$self->input_files};
  die("There can only ever be one bam file associated with a result set!") unless(@one_bam == 1);
  my $bam_file = $one_bam[0];

  my $cmd = qq(samtools index $bam_file);
  run_system_cmd($cmd, undef, 1);

  my $count_cmd = 'samtools idxstats '.$bam_file.' | awk \'{total = total + $2} END{print total}\'';
  $self->helper->debug(2, "Running:\n\t".$count_cmd);
  my $total_mapped = run_backtick_cmd($count_cmd);

  if(! $total_mapped) {
    $self->throw_no_retry("Failed to get number of mapped reads from index of:\n\t".$bam_file);
  }
  my $wiggle_tools_cmd = 'mean scale '.(10**9 / $total_mapped).' '.$bam_file;

  return $wiggle_tools_cmd;
}


1;
