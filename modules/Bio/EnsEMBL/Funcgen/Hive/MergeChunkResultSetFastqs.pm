=pod

=head1 NAME

Bio::EnsEMBL::Hive::Funcgen::MergeChunkResultSetFastqs

=head1 DESCRIPTION

This analysis validates the composition of the ResultSet and merges
and chunks control or signal fastq files appropriately, in preparation
for running the alignment of individual chunks.

=cut

package Bio::EnsEMBL::Funcgen::Hive::MergeChunkResultSetFastqs;

use warnings;
use strict;

use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( is_gzipped run_system_cmd );

use base ('Bio::EnsEMBL::Funcgen::Hive::BaseDB');

#TODO... use and update the tracking database dependant on no_tracking...

#todo
# Add status tracking support for identifying which InputSubsets might already have been aligned?
# Is this actually relevant? The files may have been moved away? Options here are to force re-run
# or move files back from archive. We should enable archive access of bams?
# Normal mode would simply try and copy the files back? There maybe parallel processes doing trying to do this
# hence checksum is necessay here. Reruns of job may pick up complete copy later on.
# We need to copy the check sum file first.
# There is still a very slight risk of a race condition when files are tested and copied.
# Hence we probably need to flock here?
# dynamic file archiving should not run on controls, leave this as a manual post pipeline step
# to minimise these risks
# Add ignore_checksums flag, this should stil create them, jsut not check them?
#



sub fetch_input {  
  my $self = shift;
  $self->SUPER::fetch_input();
  my $rset         = $self->fetch_Set_input('ResultSet');


  #RunAligner need not know anything about the ResultSet
  #So we can validate all that here and simply pass the files through to align
  #Need to leave module validation to RunAligner
  #but can pass the params directly, rather than in ResultSet
  

  
  #get/validate merge and run_controls
  #what about dataflow onwards, this will differ
  #dependant the next analysis 
  #e.g. DefineMergedDataSet, Run_BWA_and_QC_merged, DefineReplicateDataSet or BWA_ReplicateFactory
  #or can we do this implicitly just by check if we have set_names/ids set?
  
  
  my $run_controls = $self->get_param_method('run_controls', 'silent');#change this to param_silent, do we need this in run?
  #This may allow unmerged controls, if we set merge to 0 in the analysis config
  my $merge        = $self->get_param_method('merge', 'silent', $run_controls);
  
  
  
  
  #we really need to define a Base Aligner class similar to the PeakCaller class
  #define a standard interface/requirements
  #location of indexes needs to be built here
  #based an analysis name
  #do we need a index_required sub?
  #Let's do all of this in BWA for now? 
  
  #Hmm this is just chunking the files! and submitting the individual jobs!
  #Hence we need more analyses:
  #1 to run the individual BWA jobs
  #2 to merge the alignments and perform QC!
  
  
  

  $self->get_param_method('fastq_chunk_size', 'silent', '16000000');#default should run in 30min-1h 
  #does this even support input sets yet?

  
  #We need to define the work dir here for the intermediate chunk/alignments files
  #output_dir here is for alignment (no need for repository)
  $self->get_ouput_work_dir_methods($self->alignment_dir($rset));
 
  $self->helper->debug(1, "Work dir is:\t".$self->work_dir.
                            "\nOutput dir is:\t".$self->output_dir);
  
  
    
  return 1;
}

sub run {
  my $self         = shift;
  my $rset         = $self->ResultSet;
  my $run_controls = $self->run_controls;
  my $merge        = $self->merge;
  
  
  #Maybe we need to handle pre-aligned ResultSets here
  #check status and file
  
  #unsafe to check individual input_subsets, as these may have been processed
  #in an incomplete manner previously
  #Although we need to be mindful that a ResultSet may have had an InputSubset added
  #hence we can't re-use an old alignment file
  #this should be handled when creating/rolling back the ResultSet
  
  
  #This status needs to be CS specific!!
  my $align_status = $self->get_coord_system_status('ALIGNED');#put this in BaseSequenceAnalysis
  
  #We have an inheritance issue here
  #BaseSequenceAnalysis isa BaseImporter?
  
  if($rset->has_status($align_status)){
    throw("Need to implement force/recover_alignment. Found $align_status ResultSet:\t".
      $rset->name."\n");
     
    #Actually, this should only be allowed if we are recovering
    #force should have rolled back the ResultSet ALIGNED status
    #although ideally this should be handled in the previous analysis
    #and flowed directly to DefineReplicate/MergedDataSet
    #Hence we should never reuse the merged fastq?
    #if we are recovering, we want the bam file (given the reps are the same)
    #if we are forcing or rolling back, then we should probably redo everything
  }
  
  
  

  my (@local_urls, @rep_numbers);
  my $throw = '';
  
 ISSET: foreach my $isset(@{$rset->get_support('input_subset')}){

    if($isset->is_control && 
       (! $run_controls)){
      next;    
    }
 
 
    my @tinfo = @{$self->tracking_adaptor->fetch_InputSubset_tracking_info($iset)};
    #This is currently returning an ARRAYREF, but will change once the input_subset model
    #is corrected
    #TODO check whether input_subsets are unique and tidy as required!
    
    if(! @tinfo){
      $throw .= "Could not find tracking info for InputSubset:\t".
        $isset->name."\n";
      next;
    }
    
    #todo define the align output file using the same code in get_alignment_file_by_InputSet?
    #sub out the file name bit 
  
   
      
    foreach my $tr(@{$self->tracking_info}){
      my $local_url = $tr->{local_url}; #Defined here until we update the schema
      my $found_path;
     
      if(! defined $local_url){
        $throw .= "Found an InputSubset without a local_url, has this been downloaded?:\t".
          $isset->name."\n";
        next ISSET;  
      }
  
      push @rep_numbers, $isset->replicate;
      push @local_urls, $found_path;      
    }
  } 
  
  
  throw($throw) if $throw;
  
  if((scalar(@local_urls) > 1) &&
     ! $merge){
    throw('ResultSet '.$rset->name.
      " has more than one InputSubset, but merge has not been specified:\n\t".
      join("\n\t", @local_urls));    
  }  

  
  my @fastqs;
  
  foreach my $fastq_file(@local_urls){   
    
    #This should unzip the file to do the checksum, as we shouldn't have checksums
    #on zipped files
    
    
    eval { $found_path = check_file($local_url, 'gz', {checksum => 1}); };
    
    if($@){
      $throw .= "$@\n";
      next ISSET;  
    }
    elsif(! defined $found_path){
      $throw .= "Could not find InputSubset local_url path, is either not downloaded, deleted or in warehouse:\t".
        $isset->name."\n";
          
      #Could try warehouse here?
    }
     
    push @fastqs, $fastq_file;
  }    
     
  
  
  #Issues here in that the controls may not be from this experiment/study.
  #Hence we can't get the study name from the ResultSet, we need to get if form the InputSubsets themselves
  
  #Also rename this MergeChunkResultSetFastqs?
  
  
     
  #Cat rep numbers so we know exactly what this is. 
  my $merged_file_prefix = $self->get_study_name_from_ResultSet($rset).
    '_'.join('_', sort(@rep_numbers)).'.fastq';
    
  #Clean away any that match the prefix, so we don't pick them up erroneously
  #after the split.  
  #Is an under score the correct thing to concat the prefix with the chunk suffix?
  run_system_cmd("rm -f ${merged_file_prefix}_*");
     
  my $cmd = 'cat '.join(' ', @fastqs).' | split -d -a 4 -l '.
    $self->fastq_chunk_size.' - '.$self->work_dir.'/'.$merged_file_prefix;
   
  run_system_cmd($cmd);

  #Now I need to know exactly what files were produced to data flow to individual alignment jobs
  $self->set_param_method('chunk_files', [@{run_backtick_cmd("ls ${merged_file_prefix}_*")}]);

  return;
}


sub write_output {  # Create the relevant jobs
  my $self = shift @_;

  # files to align
  $self->dataflow_output_id($self->chunk_files, 2);

  # merge data across replicates
  $self->dataflow_output_id(, 3);#input_id
  return 1;

}



1;
