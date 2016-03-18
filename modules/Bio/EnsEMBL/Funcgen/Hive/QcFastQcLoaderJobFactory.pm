
=pod 

=head1 NAME

Bio::EnsEMBL::Funcgen::Hive::QcFastQcLoaderJobFactory

=head1 DESCRIPTION

standaloneJob.pl Bio::EnsEMBL::Funcgen::Hive::QcFastQcLoaderJobFactory -input_id '{chromosome_file => "/lustre/scratch109/ensembl/funcgen/mn1/ersa/faang/reference_files/CCAT/homo_sapiens_.CCAT_chr_lengths.txt", batch_param_names => ["no_write","rollback","full_delete","slices","skip_slices","result_set_only","result_set_mode","recover","alignment_analysis","peak_analysis","permissive_peaks","control_feature","no_idr","indexed_ref_fasta","idr_analysis","max_peaks","checksum_optional"], use_tracking_db => 1, dnadb => {"-dnadb_host" => "ens-livemirror","-dnadb_name" => "homo_sapiens_core_82_38","-dnadb_pass" => "","-dnadb_port" => 3306,"-dnadb_user" => "ensro"}, out_db => {"-dbname" => "mn1_faang2_tracking_homo_sapiens_funcgen_81_38","-host" => "ens-genomics1","-pass" => "ensembl","-port" => 3306,"-user" => "ensadmin"}, pipeline_name => "blah", data_root_dir => "/lustre/scratch109/ensembl/funcgen/mn1/ersa/faang/", "alignment_analysis" => "bwa_samse","checksum_optional" => 0,"dbID" => "2270","filter_from_format" => undef,"max_peaks" => 24914,"peak_analysis" => "SWEmbl_R0005","permissive_peaks" => "SWEmbl_R0005","set_name" => "F36P:hist:BR1_H3K4me3_3526_bwa_samse","set_type" => "ResultSet"}'

=cut

package Bio::EnsEMBL::Funcgen::Hive::QcFastQcLoaderJobFactory;

use warnings;
use strict;

use base qw( Bio::EnsEMBL::Funcgen::Hive::BaseDB );

sub run {
  my $self = shift;
  
  my $tempdir = $self->param_required('tempdir');
  
  my $fastqc_summary_file = `find $tempdir -maxdepth 2 -mindepth 2 -type f -name summary.txt`;
  chomp($fastqc_summary_file);
  die("Can't find file $fastqc_summary_file!") unless(-e $fastqc_summary_file);
  
  use Bio::EnsEMBL::Hive::Utils;
  Bio::EnsEMBL::Hive::Utils->import(qw/stringify destringify/);
  
  my $input_id = destringify($self->input_id);
  
  $input_id->{fastqc_summary_file} = $fastqc_summary_file;
  
#   use Data::Dumper;
#   print Dumper($fastqc_summary_file, $input_id);
#   die();
  
  $self->dataflow_output_id($input_id, 2);
  return;
}

1;


