=pod 
=head1 NAME

Bio::EnsEMBL::Funcgen::Hive::QcFastQcLoaderJobFactory

=head1 DESCRIPTION
=cut

package Bio::EnsEMBL::Funcgen::Hive::QcFastQcLoaderJobFactory;

use warnings;
use strict;

use base qw( Bio::EnsEMBL::Funcgen::Hive::BaseDB );

use Bio::EnsEMBL::Hive::Utils;
Bio::EnsEMBL::Hive::Utils->import(qw/stringify destringify/);

sub run {
  my $self = shift;
  
  my $tempdir = $self->param_required('tempdir');
  
  my $fastqc_summary_files = `find $tempdir -maxdepth 2 -mindepth 2 -type f -name summary.txt`;
  chomp($fastqc_summary_files);
  
  if ($fastqc_summary_files eq '') {
    die("Couldn't find summary file in ${tempdir}!");
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
