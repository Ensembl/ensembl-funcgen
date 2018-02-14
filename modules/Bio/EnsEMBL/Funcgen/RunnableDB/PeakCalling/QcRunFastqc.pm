=pod 
=head1 NAME

Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::QcRunFastqc

=head1 DESCRIPTION
=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::QcRunFastqc;

use warnings;
use strict;

use base 'Bio::EnsEMBL::Hive::Process';

use Bio::EnsEMBL::Hive::Utils;
Bio::EnsEMBL::Hive::Utils->import(qw/stringify destringify/);

sub run {
  my $self = shift;
  
  my $tempdir   = $self->param_required('fastqc_tempdir');
  my $read_file = $self->param_required('read_file');
  
  my $cmd;
  my $has_failed;
  my $run_options = {
    use_bash_pipefail => 1
  };

  $cmd = qq(mkdir -p $tempdir);
  $has_failed = $self->run_system_command($cmd, $run_options);
  if ($has_failed) {
    $self->throw("The following command failed:\n" . $cmd)
  }

  $cmd = qq(timeout 2h fastqc --extract -o $tempdir $read_file);
  
  $has_failed = $self->run_system_command($cmd, $run_options);
  if ($has_failed) {
    $self->throw("The following command failed:\n" . $cmd)
  }
  
  my $fastqc_summary_files = `find $tempdir -maxdepth 2 -mindepth 2 -type f -name summary.txt`;
  chomp($fastqc_summary_files);
  
  if ($fastqc_summary_files eq '') {
    $self->throw(
      "Couldn't find summary file in ${tempdir}!\n"
      . "This likely means that fastqc failed to run on the read file: ${read_file}."
    );
  }
  
  my @fastqc_summary_file = split "\n", $fastqc_summary_files;
  
  my $input_id = destringify($self->input_id);
  
  foreach my $current_fastqc_summary_file (@fastqc_summary_file) {
  
    if (! -e $current_fastqc_summary_file) {
      die("Can't find file ${current_fastqc_summary_file}!");
    }
    $input_id->{fastqc_summary_file} = $current_fastqc_summary_file;
    $self->dataflow_output_id($input_id, 2);
  }
  return;
}

1;
