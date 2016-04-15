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
  
  my $peak_caller = $self->param('peak_caller');
  
  my $file_extension;
  
  if ($peak_caller eq 'CCAT') {
    $file_extension = 'significant.region';
  }
  if ($peak_caller eq 'SWEmbl') {
    $file_extension = 'txt';
  }
  
  die unless ($file_extension);
  
  my $peaks_file = join '/', (
    $dirname,
    $file_prefix . '.' . $file_extension
  );

  die unless(-e $peaks_file);
  
  my $input_id = $self->create_input_id(
    $feature_set, 
    $analysis, 
    $result_set,
    $peaks_file
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

  my $feature_set_id = $feature_set->dbID;
  
  my $work_dir = $self->param_required('work_root_dir');
  
  my $epigenome_production_name = $feature_set->epigenome->production_name;
  
  my $temp_dir = "$work_dir/quality_checks/ProportionOfReadsInPeaks/$epigenome_production_name/$feature_set_id";
  my $out_db = $self->param('out_db');
  my $peak_caller = $self->param('peak_caller');

=head1

 {
 "aligner_param_methods" => ["sam_ref_fai"],
 "alignment_analysis" => "bwa_samse",
 "assembly" => "GRCh38",
 "batch_param_names" => ["no_write","rollback","full_delete","slices","skip_slices","result_set_only","result_set_mode","recover","alignment_analysis","peak_analysis","permissive_peaks","control_feature","no_idr","indexed_ref_fasta","max_peaks","checksum_optional"],
 "broad_peak_feature_types" => ["H3K36me3","H3K27me3","H2AK5ac","H2BK12ac","H3K14ac","H3K23me2","H3K4me1","H3K79me1","H3K79me2","H3K9me1","H3K9me3","H4K20me1","H4K8ac"],
 "can_DefineMergedDataSet" => 1,
 "can_PreprocessIDR" => 1,
 "can_run_SWEmbl_R0005_replicate" => 1,
 "checksum_optional" => 0,
 "control_feature_types" => ["Goat-IgG","Rabbit-IgG","WCE","rat-IgG-control","rabbit-IgG-control","mouse-IgG-control","GFP"],
 "data_root_dir" => "/lustre/scratch109/ensembl/funcgen/mn1/ersa/mn1_refactorMerge2_tracking_homo_sapiens_funcgen_81_38",
 "dbID" => 1016,
 "default_gender" => "male",
 "dnadb" => {"-dnadb_host" => "ens-livemirror","-dnadb_name" => "homo_sapiens_core_83_38","-dnadb_pass" => undef,"-dnadb_port" => 3306,"-dnadb_user" => "ensro"},
 "fastq_chunk_size" => 1000000,
 "fastq_root_dir" => "/lustre/scratch109/ensembl/funcgen/mn1/ersa/mn1_refactorMerge2_tracking_homo_sapiens_funcgen_81_38/fastq",
 "feature_set_analysis_logic_name" => "ccat_histone",
 "filter_from_format" => undef,
 "hive_output_dir" => "/lustre/scratch109/ensembl/funcgen/mn1/ersa/mn1_refactorMerge2_tracking_homo_sapiens_funcgen_81_38/output/ersa4ever/hive_debug",
 "out_db" => {"-dbname" => "mn1_refactorMerge2_tracking_homo_sapiens_funcgen_81_38","-driver" => "mysql","-host" => "ens-genomics1","-pass" => "ensembl","-port" => 3306,"-user" => "ensadmin"},
 "out_db_url" => "mysql://ensadmin:ensembl\@ens-genomics1:3306/mn1_refactorMerge2_tracking_homo_sapiens_funcgen_81_38",
 "peak_caller" => "CCAT",
 "pipeline_name" => "ersa4ever",
 "set_name" => "UT7:hist:BR1_H3K27me3_3526_ccat_histone",
 "set_type" => "DataSet",
 "source" => "CallBroadPeaks",
 "species" => "homo_sapiens",
 "use_tracking_db" => 1,
 "work_root_dir" => "/lustre/scratch109/ensembl/funcgen/mn1/ersa/mn1_refactorMerge2_tracking_homo_sapiens_funcgen_81_38/output/ersa4ever"
 }
 
/lustre/scratch109/ensembl/funcgen/mn1/ersa/mn1_refactorMerge2_tracking_homo_sapiens_funcgen_81_38/output/mn1_refactorMerge2_tracking_homo_sapiens_funcgen_81_38/result_feature/

  cmd => qq( proportion_of_reads_in_peaks.pl )
  . qq( --peak_file #peak_file#              )
  . qq( --temp_dir #temp_dir#                )
  . qq( --bam_file #bam_file#                )
  . qq( --peak_caller #peak_caller#          )
  . qq( --feature_set_id #feature_set_id#    )
  . qq( --user   #tracking_db_user#   )
  . qq( --pass   #tracking_db_pass#   )
  . qq( --host   #tracking_db_host#   )
  . qq( --dbname #tracking_db_name#   )

=cut

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
