=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.


=head1 NAME

Bio::EnsEMBL::Hive::Funcgen::PreprocessFastqs

=head1 DESCRIPTION

This analysis validates the composition of the ResultSet and merges
and chunks control or signal fastq files appropriately, in preparation
for running the alignment of individual chunks.

=cut

package Bio::EnsEMBL::Funcgen::Hive::PreprocessFastqs;

use warnings;
use strict;

use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( is_gzipped run_system_cmd
                                               get_set_prefix_from_Set
                                               run_backtick_cmd check_file );

use base qw( Bio::EnsEMBL::Funcgen::Hive::BaseDB );

sub fetch_input {  
  my $self = shift;
  
  $self->SUPER::fetch_input();
  my $result_set = $self->fetch_Set_input('ResultSet');

  my $run_controls = $self->get_param_method('result_set_groups', 'silent') ? 1 : 0;
  $self->set_param_method('run_controls', $run_controls);
  my $merge        = $self->get_param_method('merge', 'silent', $run_controls); 
  $self->get_param_method('checksum_optional', 'silent');

  $self->get_param_method('fastq_chunk_size', 'silent', '16000000');#default should run in 30min-1h 
  
  $self->get_output_work_dir_methods($self->alignment_dir($result_set, 1, $run_controls));#default output_dir 
  return;
}

sub run {
  my $self         = shift;
  my $result_set   = $self->ResultSet;
  my $run_controls = $self->run_controls;
  my $merge        = $self->merge;
 
  if($run_controls) {
    # control flag
    my $exp = $result_set->experiment(1);
  }
#   use Data::Dumper;
#   $Data::Dumper::Maxdepth = 3;
#   print Dumper($result_set);
#   die;

  my @fastqs;
  my $throw = '';
  my $set_rep_suffix = '';
  my @input_subsets = @{$result_set->get_support};
  
  # The code here is trying to figure out where in the pipeline it is being 
  # run. If it realises it is splitting the fastq files in the idr step,
  # it does something special.
  #
  my $is_run_during_idr_step
    = 
      (! $run_controls) && 
      (! $merge) &&
      ($self->is_idr_FeatureType($result_set->feature_type));
  
  if ($is_run_during_idr_step) {
  
    # Get the input_subset object with the signal
    my @signal_input_subsets = grep { ! $_->is_control } @input_subsets;
    
    my %temp =  map {
      ( $_->biological_replicate => 1 )
    } @signal_input_subsets;
    my @biological_replicate_number = keys %temp;
    
    %temp =  map {
      ( $_->technical_replicate => 1 )
    } @signal_input_subsets;
    my @technical_replicate_number = keys %temp;
    
#     use Data::Dumper;
#     print Dumper(\@biological_replicate_number);
#     print Dumper(\@technical_replicate_number);
    
    # Assert there is only one.
    if(scalar(@biological_replicate_number) != 1) {
      $self->throw_no_retry('Expected 1 InputSubset(replicate) for IDR ResultSet '.
        $result_set->name.' but found '.scalar(@signal_input_subsets).' Specify merge, or restrict to 1 InputSubset');  
    }

    # This is used for creating subdirectories. Then it is passed on as a 
    # batch parameter, but probably never used in any of the following 
    # analyses.
    #
    $set_rep_suffix = 
      '_BR_' . $biological_replicate_number[0] . 
      '_TR_' . join '_', @technical_replicate_number;
#     die ($set_rep_suffix);
  }
  
  foreach my $current_input_subset (@input_subsets) {

    if(
      (  $current_input_subset->is_control && ! $run_controls) || 
      (! $current_input_subset->is_control &&   $run_controls)
    ){
      next;
    }
 
    if(! $self->tracking_adaptor->fetch_tracking_info($current_input_subset)){
       $throw .= "Could not find tracking info for InputSubset:\t".
        $current_input_subset->name."\n";
      next;
    }
    
    if(! defined $current_input_subset->local_url){
      $throw .= "Found an InputSubset without a local_url, has this been downloaded?:\t".
        $current_input_subset->name."\n";
      next;
    }

    my $found_path;
    my $params = {};

    if(defined $current_input_subset->md5sum || ! $self->checksum_optional ){
#       $params->{checksum} = $current_input_subset->md5sum;
    }
    
    my $local_url = $current_input_subset->local_url;
    #Look for gz files too. These would normally already be gzipped
    #if downloaded from a repository
    #But they may have been gzipped after processing if produced locally
    #add .tgz support here?
    #we can't do a md5 check if we don't match the url exactly
    eval { $found_path = check_file($local_url, 'gz', $params); };
 
    if($@){
      $throw .= "$@\n";
      next;  
    }
#     elsif(! defined $found_path){
#       $throw .= "Could not find fastq file, is either not downloaded, has been deleted or is in warehouse:\t".
#         $local_url."\n";
#       #Could try warehouse here?
#     }
#     elsif($found_path !~ /\.(?:t){0,1}gz$/o){
#       #use is_compressed here?
#       #This will also modify the original file! And potentially invalidate any checksumming
#       $self->throw_no_retry("Found unzipped path, aborting as gzipping will invalidate any further md5 checking:\t$found_path");
#       #run_system_cmd("gzip $found_path");
#       #$found_path .= '.gz';  
#     }
    push @fastqs, $found_path;  
  }
 
  throw($throw) if $throw;
  
#   if((scalar(@fastqs) > 1) &&
#      ! $merge){
#     throw('ResultSet '.$result_set->name.
#       " has more than one InputSubset, but merge has not been specified:\n\t".
#       join("\n\t", @fastqs));    
#   }  
 
  my $set_prefix = get_set_prefix_from_Set($result_set, $run_controls).
    '_'.$result_set->analysis->logic_name.$set_rep_suffix; 
    
  my $current_working_directory = $self->work_dir . '/' . $result_set->dbID;
  
  run_system_cmd("mkdir -p $current_working_directory");

  # For safety, clean away any that match the prefix
  # No exit flag, in case rm fails due to no old files
  #
  run_system_cmd('rm -f '.$current_working_directory."/${set_prefix}.fastq_*", 1);

  my $cmd = 'zcat ' . join(' ', @fastqs) . ' | perl -pe \'s/\@SRR[^ ]+? /@/g\' | ' . ' split --verbose -d -a 4 -l ' 
    . $self->fastq_chunk_size . ' - ' . $current_working_directory . '/' . $set_prefix.'.fastq_';

  $self->helper->debug(1, "Running chunk command:\n$cmd");
  
  my @split_stdout = run_backtick_cmd($cmd);
  (my $final_file = $split_stdout[-1]) =~ s/creating file \`(.*)\'/$1/;
  
  if(! defined $final_file) {
    $self->throw_no_retry('Failed to parse (s/.*\`([0-9]+)\\\'/$1/) final file '.
      ' from last split output line: '.$split_stdout[-1]);  
  }
  
  # Get files to data flow to individual alignment jobs
  my @fastq_files = run_backtick_cmd('ls '.$current_working_directory."/${set_prefix}.fastq_*");
  @fastq_files    = sort {$a cmp $b} @fastq_files;
  
  # Now do some sanity checking to make sure we have all the files
  if($fastq_files[-1] ne $final_file) {
    $self->throw_no_retry("split output specified last chunk file was numbered \'$final_file\',".
      " but found:\n".$fastq_files[-1]);  
  } else {
    $final_file =~ s/.*_([0-9]+)$/$1/;
    $final_file =~ s/^[0]+//;
    $final_file ||= 0;  #Handle the 0000 case
  
    $self->debug(1, "Matching final_file index $final_file vs new_fastq max index ".$#fastq_files);
    
    if($final_file != $#fastq_files){
      $self->throw_no_retry('split output specified '.($final_file+1).
        ' file(s) were created but only found '.scalar(@fastq_files).":\n".join("\n", @fastq_files));  
    }
  }
  
  $self->set_param_method('fastq_files', \@fastq_files);
  my %batch_params = %{$self->batch_params};
  
  # Create funnel job first so hive can start creating the jobs in the database.
  #
  my %signal_info;  
  if($run_controls) {
    # for flow to MergeControlAlignments_and_QC
    %signal_info = (result_set_groups => $self->result_set_groups);
  }
 
  foreach my $fastq_file(@{$self->fastq_files}) {
    $self->dataflow_output_id(
      {
	%batch_params,
	# We could regenerate this from result_set and run controls
	output_dir => $current_working_directory, 
	gender     => $result_set->epigenome->gender,
	analysis   => $result_set->analysis->logic_name,   
	query_file => $fastq_file
      }, 2
    );
  }
  
  my $file_prefix  = $self->get_alignment_path_prefix_by_ResultSet($result_set, $run_controls); 
  
  my $bam_file_with_unmapped_reads_and_duplicates = $file_prefix.'.with_unmapped_reads_and_duplicates.bam';
  
  # Name of the final deduplicated bam file
  #
  my $bam_file = $file_prefix . '.bam';

  $self->dataflow_output_id(
    {
	%batch_params,
	set_type   => 'ResultSet',
	set_name   => $result_set->name,
	dbID       => $result_set->dbID,
	fastq_files => $self->fastq_files,
	output_dir => $self->output_dir,
	set_prefix => $set_prefix,
	bam_file_with_unmapped_reads_and_duplicates => $bam_file_with_unmapped_reads_and_duplicates,
	bam_file => $bam_file,
	%signal_info
      }, 3
  );

  return;
}

1;
