package Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::RunPermissiveSWEmbl;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;

sub run {

  my $self = shift;
  
  my $data_root_dir = $self->param_required('data_root_dir');
  my $plan          = $self->param_required('plan');
  
  my $bam_file_no_duplicates = $plan
    ->{remove_duplicates}
    ->{bam_file}
    ->{real}
  ;
  
  my $swembl_parameters = "-f 150 -m 8 -p 0.04 -P 0.5 -d 70 ";
  my $cmd = "SWEmbl -F $swembl_parameters -i $data_root_dir/$bam_file_no_duplicates -o $data_root_dir/${bam_file_no_duplicates}.swembl.txt";
  my $error_occurred = $self->run_system_command($cmd);
  
  if ($error_occurred) {
    $self->throw("There was an error running:\n$cmd");
  }
}

1;
