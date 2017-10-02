package Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::WriteBigWig;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;

sub run {

  my $self                = shift;
  my $species             = $self->param_required('species');
  my $execution_plan_list = $self->param_required('execution_plan_list');
  my $data_root_dir       = $self->param_required('data_root_dir');
  
  my $bam_file = $execution_plan_list
    ->[0]
    ->{call_peaks}
    ->{control_alignment}
    ->{source}
    ->{remove_duplicates}
    ->{bam_file}
    ->{real};

  my $bigwig_file = $execution_plan_list
    ->[0]
    ->{call_peaks}
    ->{control_alignment}
    ->{source}
    ->{remove_duplicates}
    ->{bigwig_file}
    ->{real};
  
  my $epigenome_gender = $execution_plan_list
    ->[0]
    ->{call_peaks}
    ->{control_alignment}
    ->{source}
    ->{remove_duplicates}
    ->{align}
    ->{to_gender};
  
  my $assembly = $execution_plan_list
    ->[0]
    ->{call_peaks}
    ->{control_alignment}
    ->{source}
    ->{remove_duplicates}
    ->{align}
    ->{to_assembly};
  
  my $bam_file_full_path    = "$data_root_dir/$bam_file";
  my $bigwig_file_full_path = "$data_root_dir/$bigwig_file";

  use Bio::EnsEMBL::Funcgen::Hive::RefBuildFileLocator;
  my $bwa_index_locator = Bio::EnsEMBL::Funcgen::Hive::RefBuildFileLocator->new;
  
  my $chromosome_lengths_relative = $bwa_index_locator->locate({
    species          => $species,
    epigenome_gender => $epigenome_gender,
    assembly         => $assembly,
    file_type        => 'chromosome_lengths_by_species_assembly',
  });
  my $reference_data_root_dir = $self->param('reference_data_root_dir');

  my $chromosome_length_file = $reference_data_root_dir . '/' . $chromosome_lengths_relative;
  
  $self->say_with_header("bam_file    = $data_root_dir/$bam_file");
  $self->say_with_header("bigwig_file = $data_root_dir/$bigwig_file");
  
  my $error_occurred;
  my $cmd;
  
  $cmd = "samtools index $bam_file_full_path";
  $error_occurred = $self->run_system_command($cmd);
  if ($error_occurred) {
    $self->throw("The following command failed:\n$cmd");
  }

  my $cmd = "samtools idxstats $bam_file_full_path | awk '{total = total + \$3} END { print total }' >&2";
  (
    $error_occurred,
    my $total_mapped
  ) = 
    $self->run_system_command(
        $cmd,
        {
            use_bash_pipefail => 1,
        }
    );
  if ($error_occurred) {
    $self->throw("The following command failed:\n$cmd");
  }
  my $wiggle_tools_cmd = 'mean scale '.(10**9 / $total_mapped).' '.$bam_file_full_path;
  
  $self->say_with_header("total_mapped = $total_mapped");
  
  use File::Basename qw( dirname basename );
  my $dirname = dirname($bigwig_file_full_path);

  use File::Path qw(make_path remove_tree);
  make_path($dirname);

  
  my $cmd_wiggletools = "wiggletools write - ". $wiggle_tools_cmd;
  my $cmd_wigToBigWig = 'wigToBigWig -fixedSummaries stdin ' . $chromosome_length_file . ' ' . $bigwig_file_full_path;
  
  my $cmd = $cmd_wiggletools . ' | ' . $cmd_wigToBigWig;
   
  $self->say_with_header($cmd);
  $error_occurred = $self->run_system_command($cmd);
  if ($error_occurred) {
    $self->throw("The following command failed:\n$cmd");
  }
}

1;
