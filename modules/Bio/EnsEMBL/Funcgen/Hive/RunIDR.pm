
=head1 LICENSE

Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

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


sub fetch_input {   # fetch parameters...
  my $self = shift;
  $self->SUPER::fetch_input;
  
  if($self->param_required('set_type') ne 'FeatureSet'){
    $self->throw_no_retry("Param set_type should be 'FeatureSet', not:\t".$self->param('set_type'));  
  }
  
  my $fset_ids = $self->get_param_method('dbIDs',  'required');
  assert_ref($rset_ids, 'ARRAY', 'FeatureSet dbIDs');
  
  #IDR analysis should probably be specified as a default/pipeline_wide and batch flowable logic_name
  #This will need to be defined in BaseSequenceANalysis as it is needed by several confs
  
  $self->set_param_method('idr_analysis', 'required');
   
  return;
}


sub run {   # Check parameters and do appropriate database/file operations... 
  my $self         = shift;
  #my $helper    = $self->helper;
  
  my $dset_adaptor = $self->out_db->get_DataSetAdaptor;
  my $idr_analysis = &scalars_to_objects($self->out_db, 'Analysis',
                                                        'fetch_by_logic_name',
                                                        [$self->idr_analysis])->[0]; 
  my $dbids = $self->dbIDs;
          
  #sanity check we have different fsets
  if($dbids->[0] == $dbids->[1]){
    $self->throw_no_retry("Pre-IDR FeatureSets are identical, dbIDs:\t".join(' ', @$dbids));  
  }                              
                       
  my $fsets = &scalars_to_objects($self->out_db, 'FeatureSet',
                                                 'fetch_by_dbID',
                                                 $dbids);
  if(scalar (@$fsets) != 2){
    $self->throw("RunIDR expect 2 replicate FeatureSets:\t".
      join(' ', map {$_->name} @$fsets)));  
  }
                                               
  my $exp_name;
  my $anal_name;
  
  foreach my $fset(@$fsets){
    
    my $rsets = $dset_adaptor->fetch_by_product_FeatureSet($fset)->get_supporting_set('result');
    
    if(scalar (@$rset) != 1){
      $self->throw_no_retry("Could not find unique supporting ResultSet for FeatureSet:\t".$fset->name);  
    }
    
    $exp_name  ||= $rsets->[0]->get_Experiment->name;
    $anal_name ||= $fset->analysis->logic_name;
    
    if($exp_name ne $rsets->[0]->get_Experiment->name){
      $self->throw_no_retry("IDR Replicate FeatureSets are not from the same experiment:\t".
        $fsets->[0]->name.' '.$fsets->[1]->name);  
    }
    
    if($anal_name ne $fset->analysis->logic_name){
      $self->throw_no_retry("IDR Replicate FeatureSets do not have the same analysis:\t".
        $fsets->[0]->name.' '.$fsets->[1]->name);  
    }
  }
   
  my $max_peaks = 300000;
   
  if($anal_name !~ /swembl/io){
    #$self->input_job->transient_error( 0 );
    $self->throw_no_retry('Pre-IDR peak filtering and reformating is currently only been optimised by the SWEmbl analysis');
    #This maybe fine, but not for MACS2, see IDR docs here 
    #https://sites.google.com/site/anshulkundaje/projects/idr#TOC-CALL-PEAKS-ON-INDIVIDUAL-REPLICATES
  }
  elsif($anal_name =~ /macs/io){
    $max_peaks = 100000;
    warn "Reseting max filtered peaks value to 100000 for MACS analysis:\t$anal_name\n";   
  } 
                                                 
  #This is also done in RunPeaks, so we really need a single method to do this?
  $self->get_output_work_dir_methods( $self->db_output_dir.'/peaks/'.$exp_name.'/'.$anal_name);
  my $out_dir = $self->output_dir;
  
  #outfile is currently defined by PeakCaller::init_files
  #based on the out_dir, outfile_prefix and the output_format
  #We should probably pass this whole path, such that we can centralise how the path is generated?
  #If we let the default PeakCaller formats be used, we can't know what the suffix will be at this point
  
  my %bed_files = 
   (
    pre_idr   => [$out_dir.'/'.$fset->[0]->name.'.bed', $out_dir.'/'.$fset->[1]->name.'.bed'],
    np_idr    => [$out_dir.'/'.$fset->[0]->name.'.np_idr.bed', $out_dir.'/'.$fset->[1]->name.'.np_idr.bed'],
    num_peaks => [],
   ); 
 
   
  #If you started with ~150 to 300K relaxed pre-IDR peaks for large genomes (human/mouse), 
  #then threshold of 0.01 or 0.02 generally works well. 
  #If you started with < 100K pre-IDR peaks for large genomes (human/mouse), 
  #then threshold of 0.05 is more appropriate. 
  #This is because the IDR sees a smaller noise component and the IDR scores get weaker. 
  #This is typically for use with peak callers that are unable to be adjusted to call large number of peaks (eg. PeakSeq or QuEST)
  #What exactly are we counting here? Total number peaks across rep, average, or the max between reps?
  #This also depends on the prevalence of the mark, it may be that a particular feature type genuinely does not have many genome wide hits
    
  my ($cmd, $max_peaks, $idr_threshold);
   
  for my $bed_idx(0, 1){
    $cmd = 'wc -l '.$bed_files{pre_idr}[$bed_idx];
    $bed_files{num_peaks}[$bed_idx] = run_backtick_cmd($cmd) - 14;
  } 
  
  if( ($bed_files{num_peaks}[0] < $max_peaks) ||
      ($bed_files{num_peaks}[1] < $max_peaks) ){
        
   
  
  
    my $idr_threshold_0 = ($bed_files{num_peaks}[0] < 100000) ? 0.05 : 0.01;
    my $idr_threshold_1 = ($bed_files{num_peaks}[1] < 100000) ? 0.05 : 0.01;
    
    if($num_peaks_1 != $num_peaks_0){
      #$self->input_job->transient_error( 0 );
      $self->throw_no_retry("Identified different optimal thresholds due to pre-IDR peak counts:\n".
        "\t".$bed_files{pre_idr}[0]."\t".$bed_files{num_peaks}[0]."(${num_peaks_0})\n".
        "\t".$bed_files{pre_idr}[1]."\t".$bed_files{num_peaks}[1]."(${num_peaks_1})\n".);
    }
    
    $idr_threshold = $idr_threshold_0;
    
    $max_peaks = ($bed_files{num_peaks}[0] < $bed_files{num_peaks}[1]) ? 
                   $bed_files{num_peaks}[0] : $bed_files{num_peaks}[1];
  }
  
  
  # Convert/dump narrow peak format bed files
 
  for my $bed_idx(0, 1){  
    
    #If we get between 100-150k, then we should warn?
    if( ($bed_files{num_peaks}[$bed_idx] < 150000) &&
        ($bed_files{num_peaks}[$bed_idx] > 100000) ){
      warn 'Pre-IDR peaks counts fall in threshold grey zone, defaulting to 0.01'.
        " but maybe consider 0.02 for:\n\t".$bed_files{num_peaks}[$bed_idx]."\n";
    }
    
    
    #NarrowPeak files are in BED6+4 format:
    #chrom, chromStart, chromEnd, name, score, strand, signalValue, pValue, qValue, peak
    #The narrowPeak format has 3 columns that can be used to rank peaks 
    #(signal.value, p.value (-log_10) and q.value (-log_10)).
    my $bed_file = $bed_files{swembl}->[$bed_idx];
    
    if(! -e $bed_file){  
      $self->throw_no_retry('Need to dump bed from DB directly to required format!');
    }
    else{
      #Never re-use output file in case it is truncated due to job failure/exit.
      
      #SWEmbl output header::
      #Input  GSE30558_GSM758361_PrEC.sorted.bam
      #Reference      GSE30558_GSM758360_LNCaP.sorted.bam
      #Sequence length        0
      #Fragment length        0
      #Background     0.000000
      #Position Background    0.036383
      #Long Background        0.181917
      #Threshold      5.000000
      #Minimum count above bg 15
      #Penalty increase       70
      #Quality cutoff 0.000000
      #Result cutoff  0.000000
      #Penalty factor 0.552834
      
      #and fields:
      #Region        - Part of the genome build e.g. chromosome
      #Start pos.    - Base in region where peak starts
      #End pos.      - Base in region where peak ends
      #Count         - Number of reads in experimental sample in peak
      #Length        - Length of peak (distance between start and end pos.)‏
      #Unique pos.   - Number of unique bases within peak at which reads begin
      #Score         - The SWEMBL score, which is basically the count of filtered thresholded reads in the peak, minus the penalties (gap distances and reference sample reads).
      #Ref. count    - Number of reads in reference sample in peak
      #Max. Coverage - Depth of reads at summit
      #Summit        - Median position of highest read coverage
      
      #signalValue field was being set to (Count - Ref. Count)/min
      #Min is not really required here for ranking and simply add the header requirement
      #would be better to simply omit and skip the commented header if present?  
      #my $cmd = 'awk \'BEGIN {OFS="\t"} NR == 9 {min=$5} NR > 14 {print $1,$2,$3,".",$7,".",($4-$8)/min,-1,-1,int($9-$1)}\' '.
      
      #Now we are just skipping header and setting signal.value to score     
      $cmd = 'awk \'BEGIN {OFS="\t"} if($1 !~ /^#/) {print $1,$2,$3,".",$7,".",$7,-1,-1,int($9-$1)}\' '.
        "$bed_file | sort -k 7nr,7nr | head -n $max_peaks > ".$bed_files{np_idr}->[$bed_idx];
      run_system_cmd($cmd);    
      #Will failures of downstream pipes be caught?   

      #Sanity check we have the file with the correct number of lines
      $cmd = "wc -l $bed_files{np_idr}->[$bed_idx]";
      my $filtered_peaks = run_backtick_cmd($cmd);
      
      if($max_peak != $filtered_peaks){
        throw("Expected $max_peaks in filtered pre-IDR bed file, but found $filtered_peaks:\n\t".
          $bed_files{np_idr}->[$bed_idx]);  
      }
    }
  }
  
  #IDR analysis
  #TODO install idrCode in /software/ensembl/funcgen
  $cmd = 'Rscript ~dz1/utils/idrCode/batch-consistency-analysis.r '.
    join(' ', @{$bed_files{np_idr}})." -1 ${out_dir}/${comb_name} 0 F signal.value";  
  #signal.value is ranking measure here i.e. SWEmbl score             
  #where is output dir?                                  
  run_system_cmd($cmd);      
  
  #Do this here rather than in post_process so we parallelise the awk.
  $cmd = "awk '$11 <= $idr_threshold {print $0}' ${out_dir}/${comb_name}-overlapped-peaks.txt | wc -l";
  my $num_peaks = run_backtick_cmd($cmd);

  #Now, do we write this as an accu entry in the hive DB, or do we want it 
  #in the tracking DB?
  #Probably the later, such that we can drop/add reps to an IDR set after we have dropped the hive DB.
  #There is currently no logic place to put this in the tracking DB!
  #We would have to add a result_set_idr_stats table to handle the multiplicity
  #This is probably a good place to store the other IDR stats too?
  #Just write to file for now until we know if/what we want in the table.
  #Do we need to be concerned if thresholds differ between combinations? 
  $cmd = "echo -e \"IDR Peaks\tIDR Threshold\tPre-IDR Peaks\n$num_peaks\t$idr_threshold\t$max_peaks\" > ${out_dir}/${comb_name}-num-peaks.txt";
  run_system_cmd($cmd);

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
  #shift->dataflow_job_groups;
  return;
}

1;