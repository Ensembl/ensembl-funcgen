=pod

=head1 NAME

Bio::EnsEMBL::Hive::RunnableDB::Funcgen::SetupAlignmentPipeline

=head1 DESCRIPTION

'SetupAlignmentPipeline' Does all the setup before the Alignment is run
Checks for existence of input files, etc...
This Runnable CAN be run multiple times in parallell!

=cut

package Bio::EnsEMBL::Funcgen::Hive::Run_QC_and_Aligner;

use warnings;
use strict;

use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( is_gzipped );

use base ('Bio::EnsEMBL::Funcgen::Hive::BaseDB');

#TODO... use and update the tracking database dependant on no_tracking...

#params

#add mandatory dataflow params in here i.e. fetch_Set_input requirements

#fastq_chunk_size
#control_file 
#alignment_file
#control_feature_type
#feature_file_format
#output_dir (redefining this will require explicit definition of align_file in PreprocessAlignment)
#ignore_embargoed
#force_embargoed
#release_month (number or abrreviation)

#Plus all the params you can define for BaseDB and Base



#Do we need to support InputSet generation here?
#or should this be handled in registration?

#todo
#Move embargo date handling to IdentifySetInputs

sub fetch_input {  
  my $self = shift $_[0];
  $self->SUPER::fetch_input();


  #This analysis is currently never associated with any Set
  #ResultSets are (or should be) linked to a 'Colelction' analysis
  my $analysis = $self->out_db->get_AnalysisAdaptor->fetch_by_logic_name
                  ($self->get_param_method('alignment_analysis', 'required'));
                  
  if(! defined $analysis){
    throw("Failed to fetch alignment_analysis with logic_name:\t".$self->alignment_analysis.
          "Please select a diffent analysis or pre-load the ".
          $self->alignment_analysis.' analysis'); 
  }


  my $align_module = $analysis->module;
  eval { require $align_module; };
  
  if($@){
    throw('Failed to require alignment_analysis module:\t$align_module"); 
  }
  
  #validate program_file isn't a path?
  my $pfile_path = (defined $self->bin_dir) ? 
    $self->bin_dir.'/'.$analysis->program_file : $analysis->program_file;

  $self->set_param_method('program_file', $pfile_path);

  #Can't create runnable yet until we have pre-processed the fastq


  $self->get_param_method('fastq_chunk_size', 'silent', '16000000');#default should run in 30min-1h 
  #does this even support input sets yet?
  my $iset         = $self->fetch_Set_input('input');
  
  
  
  $self->get_ouput_work_dir_methods($self->alignment_dir($iset));
   
  #my $fq_dir = $self->validate_dir_param
  #              ('fastq_dir', 
  #               undef, #don't create 
  #               $self->param_required('fastq_root_dir').'/'.$self->species.'/'. #default
  #                 $self->get_study_name_from_InputSet($iset)
  #              ); 
  #Don't need this now we are using local_url
  
  
  #Integrating old base Alignment.pm stuff here

  #fastq input_dir are handled by local_url
  
  #We need to define the work dir here for the intermediate chunk files
  $self->get_output_work_dir_methods($self->alignment_dir);
  #output_dir here is for alignment (no need for repository)
  #we can use the work_dir here for both fastq and alignment intermediates
    
  $self->helper->debug(1, "Work dir is:\t".$self->work_dir.
                            "\nOutput dir is:\t".$self->output_dir);
  
  
  
    
  #todo define the align output file using the same code in get_alignment_file_by_InputSet?
  #sub out the file name bit 
  
  #now validate each of the InputSubSet files exist
  #if there is more than one replicate, 
 
  
  
  #todo validate we only have one here for now as we are dealing with merged sets
  #can we fix this as we are running?
  #Do we even need to throw here?
  #Can't we just handle this now?
  #as this is the step that we are actually doing the merging?
 
  #my @issets = @{$iset->get_InputSubsets};
  
  #if(scalar(@issets) != 1){
  #  throw($iset->name.' does not have 1 merged InputSubset, found '.scalar(@issets));
  #} 
 
 
  #How are going to handle control files?
  #We need to be able to track the status of individual files
  #but this is not possible at present due to the way InputSubsets are merged
  
  #How can we get around this for now without coding too many hacks?
  
  #There maybe clashes between InputSet runs, trying to align the same control file
  
  #Can we simply test for output file
 
 
  #we need this as the trackign_info doesn't contain any identifiers
  #foreach my $isset(@{$iset->get_InputSubsets}){
    #my (undef, $rep, $url, $downloaded, $avail_date, $md5sum, $is_control,
    #    $not_pooled) = @{$tracking_adaptor->fetch_InputSubset_tracking_info($isset)};
    
   #todo add local_uri in here 
     
 
 
  #Don't need this outerloop as fetch_InputSubset_tracking_info
  #also handles InputSets
  
  $self->set_param_method('tracking_info', 
                          $self->fetch_InputSubset_tracking_info($iset), 
                          'required');
      
  my (@replicates, @controls);    
 
  my $control_gzipped = 0;
  my $control_cnt     = 0
      
  foreach my $tr(@{$self->tracking_info}){
    my $local_url = $tr->{local_url}; #Defined here until we update the schema
     
    if(! defined $local_url){
      throw("Found an InputSubset without a local_url, input_subset_id = ".
             $isset->{input_subset_id});
    }
    elsif(! -f $local_url){
      throw("Found a downloaded InputSubset where the local file does not exist:\n\t".
            $local_url);  
    }
         
     #todo handle download here
     #method should probably go in tracking adaptor, or TrackingUtils?
     #download may have to be moved to a different analysis
     #as this will require a diffent resource i.e. long or basement queue!
     #This will require some branch_config
     #based on the downloaded state
     #no way of doing a nested method call!
     #and we already have some branch config in here for the aligner!
     
     #Now populate replicates and controls
     
     #todo md5 here
     
     
     #We need to be able to know whether these are merged or not
     #As the naming of the chunk files will be identical if ther are merged or not
     #This will only matter if we are running individual replicates and
     #merged Input set concurrently
     #this should never happen?
     #let's cat the rep numbers in the end of the file names for safety
     
     
     
    #Fettling with control files here could cause clashes
    #we need to try and detect the control alignment file first 
    #and create the md5 file, to mark it as being created 
    
    #do we need another analysis to run all the control files first
    #in a non-redundant fashion?
    #then we can simply throw here if it is not present and no md5?
    
    #This doesn't fix the renaming issues
    #especially as control files are currently reundant
    #and are represented in numerous input_subsets
    
    #The correct way to do this with a corrected data model
    #would be to name them as follows:
    #1 By the original name if there is only 1 rep
    #or
    #2 By the study name followed by the input_subset_ids 
    #(no need for replicate)
    
    
    #is it time to do the patch?
    #at least for control files
    
    
    
    
    
    
    
    
     
    if($tr->{is_control}){
      
      #Need to try and detect aligned file here
      #or status
      #These will usually be shared between experiments of a study
      #Let's rename it to the study name if there is more than one replicate?
      #This may cause clashes with out controls within the same study
      #can we suffix them with input_subset_id?
      
      #This will currently clash with replicate naming
      #and subsets are merged at present!!
      #
      
      if(is_gzipped($local_url)){
        $control  
      }
      else{
        
      }
      
      push @controls, [$local_url, $tr->{replicate}];
    }
    else{
      push @replicates, [$local_url, $tr->{replicate}];  
    }    
  }
    

  $self->set_param_method('control_files', \@control_files);
  $self->set_param_method('replicates',    \@replicates);


  

  #Now let's get input for new stand along RunBWA.pm

  #bwa_index and bin needs to be passed to the runnable?
  #This means this wrapper will not be generic!

 #my $bwa_index = $self->_work_dir()."/bwa_indexes/".$species."/".$species."_".$gender."_".$assembly."_unmasked.fasta";

  #$self->_bwa_index($bwa_index);
#  my $bwa_bin = $self->param('bwa_bin') || $self->_bin_dir().'/bwa';
 # $self->_bwa_bin($bwa_bin);
  return 1;
}

sub run {
  my $self = shift @_;

  my @output_ids;
  my $set_name = $self->input_set->name;
  
  #Issues around naming the control file, as this will be shared across multiple
  #input sets
  #we need to check for the expected output file and md5
  #This will still not prevent competing jobs?
  #we can touch the md5 file first, if it is present and empty a job is running
  #else, if it is present and have an md5, then we don't need to run that alignment 
  
  #Leave out bzip support now as we don't use it or fully support it.

  foreach my $file_info (@{$self->_input_files()}){
    
    my $file       = $file_info->{'path'};
    my $replicate  = $file_info->{'replicate'};
    my $file_index = $file_info->{'file_index'};

    my $cmd;

    if($file =~ /^(.*.fastq).gz$/){
      $cmd = "gunzip -c";
    }
    elsif($file =~ /^(.*.fastq).bz2$/){
      $cmd = "bunzip2 -c"
    }
    else {
      $cmd = "cat";
    }

    $cmd .= ' '.$file.' | split -d -a 4 -l '.$self->fastq_chunk_size.' - '. 
      $self->_output_dir().'/'.$set_name."_".$replicate.'_'.$file_index.'_';

    if(system($cmd) != 0){ throw "Problems running $cmd";  }
  }


  return 1;
}


sub write_output {  # Create the relevant job
  my $self = shift @_;

  

  #OLD CODE#

  my $set_name = $self->_set_name;

  my (@align_output_ids, @merge_output_ids);

  opendir(DIR,$self->_output_dir());
  for my $split_file ( grep(/^${set_name}_\d+_\d+_\d+$/,readdir(DIR)) ){
  my $output = eval($self->input_id);
  $output->{input_file} = $split_file;
  push @align_output_ids, $output;
  }
  closedir(DIR);

  # merge data for each replicate

  for my $rep (@{$self->_replicates}){
  my $output = eval($self->input_id);
  $output->{replicate} = $rep;
  push @merge_output_ids, $output;
  }


  # files to align
  $self->dataflow_output_id(\@align_output_ids, 1);

  # merge data across replicates
  $self->dataflow_output_id($self->input_id, 2);#input_id
  return 1;

}

#Private getter / setter to the fastq chunk size
sub _fastq_chunk_size {
  return $_[0]->_getter_setter('fastq_chunk_size',$_[1]);
}

#Private getter / setter to the output_ids list
sub _output_ids {
  return $_[0]->_getter_setter('output_ids',$_[1]);
}

#Private getter / setter to the output_ids list
sub _replicates {
  return $_[0]->_getter_setter('replicates',$_[1]);
}

sub _input_files {
  return $_[0]->_getter_setter('input_files',$_[1]);
}

1;
