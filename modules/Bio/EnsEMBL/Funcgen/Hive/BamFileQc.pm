package Bio::EnsEMBL::Funcgen::Hive::BamFileQc;

use warnings;
use strict;

use base qw( Bio::EnsEMBL::Funcgen::Hive::BaseDB );

sub run {
  my $self = shift;
  
  my $result_set;
  my $set_type = $self->param_required('set_type');
  if ($set_type eq 'ResultSet') {
    $result_set = $self->fetch_Set_input('ResultSet'); 
  } else {
    my $fset     = $self->fetch_Set_input('FeatureSet');
    my $analysis = $fset->analysis;
    $result_set  = $self->ResultSet; 
  }

  my $result_set_id = $result_set->dbID;
  
  my $is_control = $self->param('is_control') ? 1 : undef;
  my $align_prefix   = $self->get_alignment_path_prefix_by_ResultSet($result_set, $is_control, 1);
  
  my $has_duplicates = $self->param('has_duplicates');
  
  my $bam_file;
  if ($has_duplicates) {
    $bam_file = $align_prefix   . '.with_duplicates.bam';
  } else {
    $bam_file = $align_prefix   . '.bam';
  }
  
  if (! -e $bam_file) {
  
    use Bio::EnsEMBL::Utils::Logger;
    my $logger = Bio::EnsEMBL::Utils::Logger->new();
    $logger->error(
      "The bam file $bam_file doesn't exist. This can happen, if it has just "
      . "been created and the file system has not been updated yet on all "
      . "nodes. In cases like these it may be resolved by restarting after "
      . "waiting a bit. This job will sleep for 10 seconds now and then die. "
      . "Hopefully upon retry it will work.\n"
    );
    sleep(10);
    die("$bam_file doesn't exist. Dying now and hoping for more luck in "
      . "the next life.");
  }

#   my $param_hash = $self->input_job->hive_pipeline->params_as_hash;

  use Bio::EnsEMBL::Hive::Utils ('stringify', 'destringify');
  my $param_hash = destringify($self->input_id);
  
  
#    use Data::Dumper;
#    die(Dumper($param_hash));
#   die($input_id);
  
  $param_hash->{bam_file} = $bam_file;
  $self->dataflow_output_id($param_hash, 2);
  
  return;
}

1;


