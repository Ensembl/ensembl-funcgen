
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

Bio::EnsEMBL::Funcgen::Hive::PreprocessIDR

=head1 DESCRIPTION

This is a very simple module which calculates the IDR threshold and number of peaks to filter before
submitting the individual IDR jobs. This enables the correct semaphore config (fan/funnel), which is currently 
impossible from the previous analyses.

# This could be omitted in multi-sempahores are supported e.g.
#    '2->A'    => ['DefineReplicateDataSet'], #fan
#    'A->3->B' => ['RunIDR'],                 #fan and immediate funnel
#    'B->4'    => ['PostProcessIDRReplicates],

=cut


#TODO
# 1 Enable theshold override
# 2 Enable dropping ResultSets (alter job id from env). Need to list rset dbids in output?
# 3 Fan filter/sort jobs, this will require a dummy analysis and maybe a job factory
# 4 Pick up peak QC failures from status table(no FeatureSet)? Or accu?
#   Or tracking table? Accu is probably enough here
#   But we do need to track the QC for individual reps, so we can mark them as dodgy   

package Bio::EnsEMBL::Funcgen::Hive::PreprocessIDR;

use warnings;
use strict;
 
use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( scalars_to_objects run_backtick_cmd run_system_cmd );
use base ('Bio::EnsEMBL::Funcgen::Hive::BaseDB');

sub fetch_input {   # fetch parameters...
  my $self = shift;
  $self->SUPER::fetch_input;
  
  if($self->param_required('set_type') ne 'ResultSet'){
    throw("Param set_type should be 'ResultSet', not:\t".$self->param('set_type'));  
  }
  
  my $rset_ids = $self->get_param_method('dbIDs',  'required');
  assert_ref($rset_ids, 'ARRAY', 'ResultSet dbIDs');
   
  $self->set_param_method('peak_analysis', &scalars_to_objects($self->out_db, 'Analysis',
                                                               'fetch_by_logic_name',
                                                               $self->param_required('permissive_peaks')));
  return;
}


sub run {   # Check parameters and do appropriate database/file operations... 
  my $self  = shift;  
  my $rsets = &scalars_to_objects($self->out_db, 'ResultSet',
                                                 'fetch_by_dbID',
                                                 $self->dbIDs);
                                                 
  if(scalar(@$rsets) <2){
    throw("Cannot run IDR with less than 2 replicate ResultSets:\t".
      join("\t", (map {$_->name} @$rsets)));  
  }                                               
  
  my $peak_analysis = $self->peak_analysis;                                                 
  my $batch_params  = $self->batch_params;
  my $exp_name      = $rsets->[0]->get_Experiment->name;
  #This is also done in RunPeaks, so we really need a single method to do this?
  $self->get_output_work_dir_methods( $self->db_output_dir.'/peaks/'.$exp_name.'/'.$peak_analysis);
  my $out_dir = $self->output_dir;
 
  # Do this here so we don't have clashes between RunIDR jobs 

  # Validate analysis
  my $max_peaks     = 300000;
  
  if($peak_analysis !~ /swembl/io){
    #$self->input_job->transient_error( 0 );
    $self->throw_no_retry('Pre-IDR peak filtering and reformating is currently only been optimised by the SWEmbl analysis');
    #This maybe fine, but not for MACS2, see IDR docs here 
    #https://sites.google.com/site/anshulkundaje/projects/idr#TOC-CALL-PEAKS-ON-INDIVIDUAL-REPLICATES
  }
  elsif($peak_analysis =~ /macs/io){#future proofing as will currently never be tested
    $max_peaks = 100000;
    warn "Reseting max filtered peaks value to 100000 for MACS analysis:\t$peak_analysis\n";   
  } 
      
  # Validate, count and define IDR threshold 
  #If you started with ~150 to 300K relaxed pre-IDR peaks for large genomes (human/mouse), 
  #then threshold of 0.01 or 0.02 generally works well. 
  #If you started with < 100K pre-IDR peaks for large genomes (human/mouse), 
  #then threshold of 0.05 is more appropriate. 
  #This is because the IDR sees a smaller noise component and the IDR scores get weaker. 
  #This is typically for use with peak callers that are unable to be adjusted to call large number of peaks (eg. PeakSeq or QuEST)
  #What exactly are we counting here? Total number peaks across rep, average, or the max between reps?
  #This also depends on the prevalence of the mark, it may be that a particular feature type genuinely does not have many genome wide hits
  my (%pre_idr_beds, $bed_file, $cmd, $num_peaks, $lt_100k, $mt_100k, 
      $log_txt, @rep_nums, $ctrl_ids);
  
  foreach my $rset(@$rsets){
 
    if($exp_name ne $rset->get_Experiment->name){
      $self->throw_no_retry("IDR Replicate ResultSets are not from the same experiment:\t".
        $rsets->[0]->name.' '.$rset->name);  
    }
    
    my $seen_rep = 0;
    
    if($rset->table_name ne 'input_subset'){
      $self->throw_no_retry('Expected input_subset support but found '.
        $rset->table_name." support for ResultSet:\t".$rset->name);  
    }
    
    my @issets = @{$rset->get_support};
    my @ctrl_ids;
    
    foreach my $isset(@issets){
      
      if($isset->is_control){
        push @ctrl_ids, $isset->dbID;;        
      }
      else{
        if($seen_rep){
          $self->throw_no_retry("Found more than 1 replicate (non-control) InputSet supporting an an IDR ResultSet:\n\t".
            join("\n\t", map {$_->name} @issets)."\n");  
        }  
        
        push @rep_nums, $isset->replicate;
        $seen_rep = 1;
      }
    }
    
    if(! defined $ctrl_ids){
      $ctrl_ids = join(' ', (sort {$a <=> $b} @ctrl_ids));  
    }
    
    my $current_ctrl_ids = join(' ', (sort {$a <=> $b} @ctrl_ids));  
    
    if($current_ctrl_ids ne $ctrl_ids){
      $self->throw_no_retry("Found mismatched controls between IDR ResultSets:\n\t".
        join("\t", (map {$_->name} @$rsets)));  
    }
    
    
    #todo validate ResultSet analysis here too?
      
    #bed file is currently defined by PeakCaller::init_files
    #based on the out_dir, outfile_prefix and the output_format
    #We should probably pass this whole path, such that we can centralise how the path is generated?
    #If we let the default PeakCaller formats be used, we can't know what the suffix will be at this point
    $bed_file = $out_dir.'/'.$rset->name.'.bed';
    
    if(! -f $bed_file){  
      $self->throw_no_retry("Could not find pre-IDR bed file:\t$bed_file\n".
        'Need to dump bed from DB directly to required format!');
    }
    
    $cmd                     = "wc -l $bed_file";
    $num_peaks               = run_backtick_cmd($cmd) - 14;
    $pre_idr_beds{$bed_file} = $num_peaks;
   
    if($num_peaks > 100000){
      $mt_100k = 1;
      
      if($num_peaks < 150000){
        warn 'Pre-IDR peaks counts fall in threshold grey zone(100k-150k), defaulting to 0.01'.
          " but maybe consider 0.02 for:\n\t$bed_file\n";  
      }  
    }
    else{
      $lt_100k = 1;
    }

    if($num_peaks < $max_peaks){
      $max_peaks = $num_peaks;  
    }
    
    $log_txt .= $bed_file."\t".$num_peaks."\n";
  }
    
  
  #Note this does not yet support MACS yet, should prbably just ignore it as we filter to 100000
  
  if($lt_100k && $mt_100k){
    $self->throw_no_retry("Identified different optimal thresholds due to pre-IDR peak counts spanning the 100k limit:\n".
      $log_txt);    
  }
  
  #TODO We need some mechanism to restart this job, to force the threshold, or by dropping 1/some of the replicates. 
  
  my $idr_threshold = ($num_peaks < 100000) ? 0.05 : 0.01;
  my $idr_name      = $exp_name.'_'.$peak_analysis.'_'.join('_', sort @rep_nums);
  $cmd = "echo -e \"Pre-IDR File\tNum Peaks\n$log_txt\nIDR Threshold = $idr_threshold\" > ${out_dir}/${idr_name}-idr_stats.txt";
  run_system_cmd($cmd);
 
  #TODO parallelise the filtering and reformating to speed things up, so let's semphaore than to a simple CMD job.
  #can we even do this as we already have a semaphore on the RunIDR and
  #maybe with a job factory? I think this is not possible without another analysis
  #but we could use a dummy? which then submit the RunIDR and semaphores the PostProcessIDR
  #Just do here for now
  my @np_bed_files;

  foreach my $bed_file(keys %pre_idr_beds){
    (my $np_bed_file = $bed_file) =~ s/\.bed$/.np_idr.bed/;  
    #Never re-use np_idr output file in case it is truncated due to job failure/exit.
    
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
    #Length        - Length of peak (distance between start and end pos.)
    #Unique pos.   - Number of unique bases within peak at which reads begin
    #Score         - The SWEMBL score, which is basically the count of filtered thresholded reads in the peak, minus the penalties (gap distances and reference sample reads).
    #Ref. count    - Number of reads in reference sample in peak
    #Max. Coverage - Depth of reads at summit
    #Summit        - Median position of highest read coverage
    
    #signalValue field was being set to (Count - Ref. Count)/min
    #Min is not really required here for ranking and simply add the header requirement
    #would be better to simply omit and skip the commented header if present?  
    #my $cmd = 'awk \'BEGIN {OFS="\t"} NR == 9 {min=$5} NR > 14 {print $1,$2,$3,".",$7,".",($4-$8)/min,-1,-1,int($9-$1)}\' '.
      
    #Now we are stripping out the header and setting signal.value to score 
    #before sorting on score and filtering based on $max_peaks  
    #and resorting based on position  
    $cmd = 'awk \'BEGIN {OFS="\t"} if($1 !~ /^#/) {print $1,$2,$3,".",$7,".",$7,-1,-1,int($9-$1)}\' '.
      "$bed_file | sort -k 7nr,7nr | head -n $max_peaks | sort -k 1,2n > ".$np_bed_file;
    run_system_cmd($cmd);    
    #Will failures of downstream pipes be caught?   

    #Sanity check we have the file with the correct number of lines
    $cmd = "wc -l $np_bed_file";
    my $filtered_peaks = run_backtick_cmd($cmd);
      
    if($max_peaks != $filtered_peaks){
      throw("Expected $max_peaks in filtered pre-IDR bed file, but found $filtered_peaks:\n\t".$np_bed_file);  
    }   
 
    
    
    #TODO check the feature_set_stat or statuses
    #such that we know the peak calling jobs has finished and passed QC!    
    #Do this for each before we submit IDR jobs, as we may have to drop some reps
    push @np_bed_files, $np_bed_file;
  }


  my $idr_analysis = &scalars_to_objects($self->out_db, 'Analysis',
                                                        'fetch_by_logic_name',
                                                        [$self->idr_analysis])->[0]; 



  #Build 2 way rep combinations for IDR jobs
  my @idr_job_ids;
  my $last_i = $#np_bed_files - 1;
  
  foreach my $i(0..$last_i){
    my $next_i = $i + 1;
  
    foreach my $j($next_i..$#np_bed_files){
      
      push @idr_job_ids, {{%{$batch_params},
                          output_dir    => $out_dir,
                          #output_prefix => $output_prefix,
                          idr_threshold => $idr_threshold,
                          idr_name      => $idr_name,
                          bed_files     => [$np_bed_files[$i], $np_bed_files[$j]]}};
    }  
  }
  
  
  
  #Now we need to pool and produce pseudo reps? This should be done way before here?!!
  #In between MergeControlAlignments_and_QC and Submit_IDR
  #
  
  
  $self->branch_job_group(2, \@idr_job_ids,
                          3, {%{$batch_params},
                              dbIDs    => $self->dbIDs,
                              set_type => 'ResultSet'});
      
  return;
}



sub write_output {  # Create the relevant jobs
  shift->dataflow_job_groups;
  return;
}

1;