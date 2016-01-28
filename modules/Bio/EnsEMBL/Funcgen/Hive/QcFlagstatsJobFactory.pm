
=pod 

=head1 NAME

Bio::EnsEMBL::Funcgen::Hive::QcFlagstatsJobFactory

=head1 DESCRIPTION

standaloneJob.pl Bio::EnsEMBL::Funcgen::Hive::QcFlagstatsJobFactory -input_id '{chromosome_file => "/lustre/scratch109/ensembl/funcgen/mn1/ersa/faang/reference_files/CCAT/homo_sapiens_.CCAT_chr_lengths.txt", batch_param_names => ["no_write","rollback","full_delete","slices","skip_slices","result_set_only","result_set_mode","recover","alignment_analysis","peak_analysis","permissive_peaks","control_feature","no_idr","indexed_ref_fasta","idr_analysis","max_peaks","checksum_optional"], use_tracking_db => 1, dnadb => {"-dnadb_host" => "ens-livemirror","-dnadb_name" => "homo_sapiens_core_82_38","-dnadb_pass" => "","-dnadb_port" => 3306,"-dnadb_user" => "ensro"}, out_db => {"-dbname" => "mn1_faang2_tracking_homo_sapiens_funcgen_81_38","-host" => "ens-genomics1","-pass" => "ensembl","-port" => 3306,"-user" => "ensadmin"}, pipeline_name => "blah", data_root_dir => "/lustre/scratch109/ensembl/funcgen/mn1/ersa/faang/", "alignment_analysis" => "bwa_samse","checksum_optional" => 0,"dbID" => "2270","filter_from_format" => undef,"max_peaks" => 24914,"peak_analysis" => "SWEmbl_R0005","permissive_peaks" => "SWEmbl_R0005","set_name" => "F36P:hist:BR1_H3K4me3_3526_bwa_samse","set_type" => "ResultSet"}'

=cut

package Bio::EnsEMBL::Funcgen::Hive::QcFlagstatsJobFactory;

use warnings;
use strict;

use base qw( Bio::EnsEMBL::Funcgen::Hive::BaseDB );

sub fetch_input {
  my $self = shift; 
  $self->SUPER::fetch_input;
  return;
}

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
  
#   use Data::Dumper;
#   $Data::Dumper::Maxdepth = 0;
#   print Dumper($input_id);
  $self->dataflow_output_id($input_id, 2);
  return;
}

sub write_output {
  my $self = shift;
  return;
}

sub create_input_id {

  my $self = shift;
  my $result_set = shift;
  my $result_set_id = $result_set->dbID;

  my $align_prefix   = $self->get_alignment_path_prefix_by_ResultSet($result_set, undef, 1);#validate aligned flag 
  my $control_prefix = $self->get_alignment_path_prefix_by_ResultSet($result_set, 1, 1);#and control flag 
  
  my $signal_bam_file  = $align_prefix   . '.bam';
  my $control_bam_file = $control_prefix . '.bam';
  
  if (! -e $signal_bam_file) {
    die("$signal_bam_file doesn't exist!");
  }
  if (! -e $control_bam_file) {
    die("$control_bam_file doesn't exist!");
  }
  
  my $out_db = $self->param('out_db');
  
  my $work_dir = $self->param_required('work_root_dir');
  
  my $temp_dir = "$work_dir/temp/flagstats/$result_set_id";
  
  use File::Path qw( make_path );
  make_path($temp_dir);
  
  use File::Basename;
  (my $signal_bam_file_base_name,  my $signal_bam_directory)  = fileparse($signal_bam_file);
  (my $control_bam_file_base_name, my $control_bam_directory) = fileparse($control_bam_file);
  
  die unless($signal_bam_directory eq $control_bam_directory);
  
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
  
  my $signal_flagstats_file  = "$temp_dir/${signal_bam_file_base_name}.signal_flagstats.txt";
  my $control_flagstats_file = "$temp_dir/${control_bam_file_base_name}.control_flagstats.txt";
  
  my $input_id = [
    {
      @$input_id_common,
      bam_file => $signal_bam_file,
      flagstats_file => $signal_flagstats_file,
      is_control => undef,
    },
    {
      @$input_id_common,
      bam_file => $control_bam_file,
      flagstats_file => $control_flagstats_file,
      is_control => 1,
    }
  ];
  
  return $input_id;
}

1;


