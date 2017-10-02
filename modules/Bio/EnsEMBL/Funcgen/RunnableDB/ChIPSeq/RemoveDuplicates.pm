package Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::RemoveDuplicates;

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
  my $bam_file_with_duplicates = $plan
    ->{remove_duplicates}
    ->{align}
    ->{bam_file}
    ->{real}
  ;
  
  eval {
    use Bio::EnsEMBL::Funcgen::Sequencing::SeqTools;
    remove_duplicates_from_bam({
      input_bam  => $data_root_dir . '/' . $bam_file_with_duplicates,
      output_bam => $data_root_dir . '/' . $bam_file_no_duplicates, 
      debug      => $self->debug,
    });
  };
  if ($@) {
    $self->throw($@);
  }
  
  #sleep(30);
  return;
}

1;
