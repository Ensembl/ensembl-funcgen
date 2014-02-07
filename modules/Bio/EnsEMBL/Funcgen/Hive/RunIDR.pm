
=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Funcgen::Hive::RunIDR

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::Hive::RunIDR;

use warnings;
use strict;
 
use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( scalars_to_objects run_system_cmd );
use base ('Bio::EnsEMBL::Funcgen::Hive::BaseDB');

### CURRENT ISSUES 
#
# 1 This is currently tied to the SWEmbl format
#
# 2 SWEmbl has not been tested for IDing large numbers of peaks, so peak filtering may have issues like MACS2
#
# 3 IDR analysis is currently not stored in the DB. Where are we going to store the threshold and num peaks used for the IDR?
#   Should probably go in the FeatureSet.description rather than different analyses, or just omit for now? 
#IDR analysis should probably be specified as a default/pipeline_wide and batch flowable logic_name
#This will need to be defined in BaseSequenceAnalysis as it is needed by several confs
#This should be merged with the SWEMBL analysis e.g. 
  

sub fetch_input {   # fetch parameters...
  my $self = shift;
  $self->SUPER::fetch_input;
 
  my $rset_ids = $self->get_param_method('bed_files',  'required');
  assert_ref($rset_ids, 'ARRAY', 'ResultSet dbIDs');
  

  $self->get_param_method('idr_threshold', 'required');
  $self->get_param_method('output_prefix', 'silent');
  $self->get_param_method('idr_name', 'silent', $self->output_prefix);
  $self->param_required('idr_name');
  $self->param_required('output_dir');
  return;
}


sub run {   # Check parameters and do appropriate database/file operations... 
  my $self          = shift;
  my $files         = $self->bed_files;
  my $output_prefix = $self->output_prefix;
  my $out_dir       = $self->output_dir;
          
  #Check we have different files
  if($files->[0] eq $files->[1]){
    $self->throw_no_retry("Pre-IDR ResultSets are identical, dbIDs:\t".join(' ', @$files));  
  }                              
 
  #Check we have 2 reps
  if(scalar (@$files) != 2){
    $self->throw_no_retry("RunIDR expect 2 replicate bed files:\t".join(' ', @$files));  
  }
       
  #Check output dir exists     
  if(! -d $out_dir){
    $self->throw_no_retry("Output directory does not exist:\t$out_dir");  
  }

  #Set output prefix if we don't have one.
  if(! defined $output_prefix){
    my @names;
    
    foreach my $file(@$files){
      (my $name = $file) =~ s/.*\///;  
      $name =~ s/(\.np_idr)*\.bed$//;
      push @names, $name;
    }
    
    $output_prefix = $names[0].'_VS_'.$names[1];      
  }

  #IDR analysis
  #TODO install idrCode in /software/ensembl/funcgen and add this an analysis?
  my $idr_name = $self->idr_name;
  
  my $cmd = 'Rscript ~dz1/utils/idrCode/batch-consistency-analysis.r '. 
    join(' ', @{$files})." -1 ${out_dir}/${idr_name} 0 F signal.value";  
  #signal.value is ranking measure here i.e. SWEmbl score                                              
  run_system_cmd($cmd);      
  
  #Do this here rather than in post_process so we parallelise the awk.
  $cmd = "awk '$11 <= ".$self->idr_threshold." {print $0}' ${out_dir}/${output_prefix}-overlapped-peaks.txt | wc -l";
  my $num_peaks = run_backtick_cmd($cmd);

  #Now, do we write this as an accu entry in the hive DB, or do we want it 
  #in the tracking DB?
  #Probably the later, such that we can drop/add reps to an IDR set after we have dropped the hive DB.
  #There is currently no logic place to put this in the tracking DB!
  #We would have to add a result_set_idr_stats table to handle the multiplicity
  #This is probably a good place to store the other IDR stats too?
  #Just write to file for now until we know if/what we want in the table.
  #Do we need to be concerned if thresholds differ between combinations? 
  
  #Warning: Parallelised appending to file!l
  $cmd = "echo -e \"IDR Comparison\tIDR Peaks\n$output_prefix\t$num_peaks\" >> ${out_dir}/${idr_name}-idr-stats.txt";
  run_system_cmd($cmd);
  
  $self->set_param_method('num_peaks', $num_peaks);
  
  
  #This probably needs writin to the accu table too! So the Postprocess jobs can pick up the min $num_peaks without reading the file
  

  #This is old 2 rep only method which sets a new threshold for a subsampled SWEmbl run 
  # as opposed to counting the number of peaks as per the IDR docs   
  #Compute cutoff value (min score) for IDR = 1%
  #my $factor = run_backtick_cmd("awk 'NF > 10 && $11 < 0.01 && (min == 0 || min > $5) {min = $5} END {print min}' ${out_dir}/${exp_name}-overlapped-peaks.txt");
  # This is when I realize that starting off in an actual scripting language would have been better... oh well
  #LOW_CUTOFF='0.0001'
  #IDR_CUTOFF='0.01'
  #new_cutoff="$(echo "$LOW_CUTOFF*$factor" | bc)"
  #my $low_cutoff = 0.0001;
  #my $new_cutoff = $factor * 
  # echo Subsample
  #samtools view -s 0.5 -b $bam1 > $sub_bam1
  #samtools view -s 0.5 -b $bam2 > $sub_bam2
  #samtools merge -f $merged_bam $sub_bam1 $sub_bam2
  #rm -f $sub_bam1 $sub_bam2
  #echo Final run with new -R param set to $new_cutoff
  #SWEMBL -F -V -i $merged_bam -f 150 -R $new_cutoff -o $output -r $control 
   
  #        $self->branch_job_group($branch, [{%{$batch_params},
  #                                           dbID       => [$rset->dbID], 
  #                                           set_name   => [$rset->name],
  #                                           set_type    => 'ResultSet'}]);){
  return;
}



sub write_output {  # Create the relevant jobs
  my $self = shift;
  #This is accumulated into an array, which is picked up by PostprocessIDR analysis on another branch.
  $self->dataflow_output_id( {'idr_peak_counts'   => $self->num_peaks}, 2);
  return;
}

1;