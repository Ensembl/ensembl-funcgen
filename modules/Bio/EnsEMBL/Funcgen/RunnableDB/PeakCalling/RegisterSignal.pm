package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::RegisterSignal;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::ExecutionPlanUtils qw (
    lock_execution_plan
);

sub run {

  my $self = shift;
  
  my $species     = $self->param_required('species');
  my $signal_plan = $self->param_required('execution_plan');
  
  lock_execution_plan($signal_plan);
  
  my $alignment_name = $signal_plan->{name};
  my $bigwig_file    = $signal_plan->{output}->{stored};

  use Bio::EnsEMBL::Funcgen::DataFile;
  my $bigwig_data_file = Bio::EnsEMBL::Funcgen::DataFile
    ->new(
      -table_id     => 0,
      -table_name   => 'alignment',
      -path         => $bigwig_file,
      -file_type    => 'BIGWIG',
      -md5sum       => undef,
    );

  my $data_file_adaptor = Bio::EnsEMBL::Registry
    ->get_adaptor(
        $species, 
        'funcgen', 
        'DataFile'
    );
  
  eval {
    $data_file_adaptor->store($bigwig_data_file);
  };
  
  if ($@) {
    my $error = $@;
    
    warn($@);
    
    my $is_already_stored_error
        = $error =~ /DBD::mysql::st execute failed: Duplicate entry/;
    
    if (!$is_already_stored_error) {
        $self->throw($@);
    }
    $bigwig_data_file = $data_file_adaptor->fetch_by_path($bigwig_file);
    
    if (! defined $bigwig_data_file) {
        $self->throw(
            "\n\nGot this when trying to store:\n\n" 
            . $@ . "\n\n"
            . "But can't fetch the data file $bigwig_file either.\n\n"
        );
    }
  }
  
  my $alignment_adaptor = Bio::EnsEMBL::Registry
    ->get_adaptor(
        $species, 
        'funcgen', 
        'Alignment'
    );
  my $alignment = $alignment_adaptor->fetch_by_name($alignment_name);
  
  $alignment->bigwig_file_id($bigwig_data_file->dbID);
  $alignment_adaptor->update($alignment);
  
  return;
}

1;
