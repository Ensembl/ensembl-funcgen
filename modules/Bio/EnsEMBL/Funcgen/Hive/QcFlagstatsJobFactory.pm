package Bio::EnsEMBL::Funcgen::Hive::QcFlagstatsJobFactory;

use warnings;
use strict;

use base qw( Bio::EnsEMBL::Funcgen::Hive::BaseDB );

sub run {
  my $self = shift;
  
  my $rset;
  my $set_type = $self->param_required('set_type');
  if ($set_type eq 'ResultSet') {
    $rset = $self->fetch_Set_input('ResultSet'); 
  } else {
    my $fset     = $self->fetch_Set_input('FeatureSet');
    my $analysis = $fset->analysis;
    $rset        = $self->ResultSet; 
  }
  
  my $input_id = $self->create_input_id($rset);
  $self->dataflow_output_id($input_id, 2);
  return;
}

sub create_input_id {

  my $self = shift;
  my $result_set = shift;
  my $result_set_id = $result_set->dbID;
  
  my $is_control = $self->param('is_control') ? 1 : undef;
  my $align_prefix   = $self->get_alignment_path_prefix_by_ResultSet($result_set, $is_control, 1);
  
  my $bam_file  = $align_prefix   . '.bam';
  
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

  my $out_db = $self->param('out_db');
  
  my $work_dir = $self->param_required('work_root_dir');
  
  my $temp_dir = "$work_dir/temp/flagstats/$result_set_id";
  
  use File::Path qw( make_path );
  make_path($temp_dir);
  
  use File::Basename;
  (my $bam_file_base_name,  my $bam_directory)  = fileparse($bam_file);
  
  my $input_id_common = [
      # Directory into which the bam files will be copied
      tempdir               => $temp_dir,
      
      result_set_id  => $result_set_id,
      
      # Connection details for the db to which the results will be written
      tracking_db_user   => $out_db->dbc->user,
      tracking_db_pass   => $out_db->dbc->password,
      tracking_db_host   => $out_db->dbc->host,
      tracking_db_name   => $out_db->dbc->dbname,
      tracking_db_port   => $out_db->dbc->port,
  ];
  
  my $flagstats_file  = "$temp_dir/${bam_file_base_name}.flagstats.txt";
  
  my $input_id = {
      @$input_id_common,
      bam_file => $bam_file,
      flagstats_file => $flagstats_file,
    };
  
  return $input_id;
}

1;


