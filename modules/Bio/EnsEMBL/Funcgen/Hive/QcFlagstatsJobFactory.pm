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
  
  my $bam_file = $self->param('bam_file_for_qc');

  my $out_db = $self->param('out_db');
  
  my $work_dir = $self->flagstats_output_dir;
  
  my $epigenome_production_name = $result_set->epigenome->production_name;
  
  my $temp_dir = "$work_dir/$epigenome_production_name/$result_set_id";
  
  use File::Path qw( make_path );
  make_path($temp_dir);
  
  use File::Basename;
  (my $bam_file_base_name,  my $bam_directory)  = fileparse($bam_file);
  
  my $flagstats_file  = "$temp_dir/${bam_file_base_name}.flagstats.txt";
  
  my $input_id = {
      bam_file       => $bam_file,
      flagstats_file => $flagstats_file,
      # Directory into which the bam files will be copied
      tempdir        => $temp_dir,
      result_set_id  => $result_set_id,
      
      # Connection details for the db to which the results will be written
      tracking_db_user   => $out_db->dbc->user,
      tracking_db_pass   => $out_db->dbc->password,
      tracking_db_host   => $out_db->dbc->host,
      tracking_db_name   => $out_db->dbc->dbname,
      tracking_db_port   => $out_db->dbc->port,
  };
  
  return $input_id;
}

1;


