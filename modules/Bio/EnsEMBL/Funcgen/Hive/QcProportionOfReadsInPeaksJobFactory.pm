package Bio::EnsEMBL::Funcgen::Hive::QcProportionOfReadsInPeaksJobFactory;

use warnings;
use strict;

use base qw( Bio::EnsEMBL::Funcgen::Hive::BaseDB );

sub run {
  my $self = shift;
  my $set_type = $self->param_required('set_type');
  
  die() unless($set_type eq 'DataSet');

  my $feature_set  = $self->fetch_Set_input('FeatureSet');
  my $analysis     = $feature_set->analysis;
  my $result_set   = $self->ResultSet; 
  
  # That is how the directory is built in RunPeaks.pm:
  # 
  # https://github.com/Ensembl/ensembl-funcgen/blob/release/83/modules/Bio/EnsEMBL/Funcgen/Hive/RunPeaks.pm#L122
  #
  my $dirname = join '/', (
    $self->peaks_output_dir,
    $result_set->experiment->name,
    $analysis->logic_name
  );
  my $file_prefix    = $result_set->name.'.'.$analysis->logic_name;
  
#   my $peak_caller = $self->param('peak_caller');
  my $peak_caller;
  my $file_extension;
  
  if ($analysis->logic_name eq 'ccat_histone') {
    $file_extension = 'significant.region';
    $peak_caller    = 'CCAT'
  } else {
    $file_extension = '.txt';
    $peak_caller    = 'SWembl'
# HACK Put better name generation into RunPeaks, then dataflow the name here.
#     $file_prefix = 'peaks';
#     $file_extension = 'swembl';
  }
  
  if (! $file_extension) {
    die;
  }
  
  my $peaks_file = join '/', (
    $dirname,
    $file_prefix . '.' . $file_extension
  );

  if (! -e $peaks_file) {
    die("Can't find peaks file $peaks_file");
  }
  
  my $input_id = $self->create_input_id(
    $feature_set, 
    $analysis, 
    $result_set,
    $peaks_file,
    $peak_caller
  );
  
  use Data::Dumper;
  print Dumper($input_id);

#   my $input_subset_id = $self->param('input_subset_id');
#   my $input_id = $self->create_input_id($input_subset_id);
#   
  $self->dataflow_output_id($input_id, 2);
  return;
}

sub create_input_id {

  my $self        = shift;
  my $feature_set = shift;
  my $analysis    = shift;
  my $result_set  = shift;
  my $peaks_file  = shift;
  my $peak_caller = shift;

  my $feature_set_id = $feature_set->dbID;
  
  my $work_dir = $self->proportion_of_reads_in_peaks_output_dir;
  
  my $epigenome_production_name = $feature_set->epigenome->production_name;
  
  my $temp_dir = "$work_dir/$epigenome_production_name/$feature_set_id";
  my $out_db = $self->param('out_db');

  # Undef means: not the control
  my $file_prefix  = $self->get_alignment_path_prefix_by_ResultSet($result_set, undef); 
  my $bam_file     = $file_prefix.'.bam';
  
  die("$bam_file does not exist!") unless(-e $bam_file);

  my $input_id = {
  
    peak_file          => $peaks_file,
    temp_dir           => $temp_dir,
    bam_file           => $bam_file,
    peak_caller        => $peak_caller,
    feature_set_id     => $feature_set_id,
    tracking_db_user   => $out_db->dbc->user,
    tracking_db_pass   => $out_db->dbc->password,
    tracking_db_host   => $out_db->dbc->host,
    tracking_db_name   => $out_db->dbc->dbname,
    tracking_db_port   => $out_db->dbc->port,

  };
  return $input_id;
}

1;
