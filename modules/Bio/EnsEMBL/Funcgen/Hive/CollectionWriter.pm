=pod 

=head1 NAME

Bio::EnsEMBL::Hive::RunnableDB::Funcgen::CollectionWriter

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::Hive::CollectionWriter;

use base ('Bio::EnsEMBL::Funcgen::Hive::BaseImporter');

use warnings;
use strict;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( generate_slices_from_names 
                                               run_system_cmd );
use Bio::EnsEMBL::Utils::Exception qw(throw);

sub fetch_input {  
  my $self = shift;

  $self->SUPER::fetch_input;    
  $self->helper->debug(1, "CollectionWriter::fetch_input after SUPER::fetch_input");  
  my $rset = $self->fetch_Set_input('ResultSet');  # Injects ResultSet, FeatureSet & DataSet methods
  $self->helper->debug(1, "CollectionWriter::fetch_input got ResultSet:\t".$rset);
  
  $self->init_branching_by_analysis;  #Set up the branch config
  my $ftype_name = $rset->feature_type->name;
 
  #$self->get_output_work_dir_methods($self->db_output_dir.'/result_feature/'.$rset->name, 1);#no work dir flag
  
  # This is required by sam_ref_fai which is called by get_alignment_file_by_ResultSet_formats
  # Move this to get_alignment_file_by_ResultSet_formats?
  $self->set_param_method('cell_type', $rset->cell_type, 'required'); 
  
  
  # TODO: need to dataflow 'bam_filtered' from alignment pipeline
  # Curretnly this is setin the alignment pipeline_wide params
  # then use this to perform the filtering in the first place
  # Of course we can data flow this, set it as a batch_param and flow it from the
  # branching analysis i.e. DefineResultSets or MergeControlAlignments_and_QC 
  # (and manually from IdentifyReplicateResultSets?)

  
  # This should be in run as it can do some conversion?
  # also need to pass formats array through 
  # dependant on which analysis we are running bam and bed for Preprocess and just bed for WriteCollections
  # Currently no way to do this dynamically as we don't know what the analysis is (?)
  # so we have to ad as analysis_params
  # and there is no way of detecting what formats down stream analyses need
  # so again, they have to be hardcoded in analysis params
  # use all formats by default
  # Can we update -input_files in the Importer after we have created it?
  
  my @file_to_delete_after_cell_line_has_been_processed;  

  if($self->FeatureSet->analysis->program eq 'CCAT') {
    
    my $exp = $rset->experiment(1);  # ctrl flag

    my $path = $self->get_alignment_path_prefix_by_ResultSet($rset, 1);
    
    my $bam_file = $path . '.bam';
    my $bed_file = $path . '.bed';
    
    if(! -e $bam_file) {
      confess("Can't find bam file $bam_file!");
    }
    
    if (! -e $bed_file) {
      my $cmd = qq(bamToBed -i $bam_file > ${bed_file}.part);
      run_system_cmd($cmd);
      $cmd = qq(mv ${bed_file}.part $bed_file);
      run_system_cmd($cmd);
      push @file_to_delete_after_cell_line_has_been_processed, $bed_file;
    }
  }
  
  my $path = $self->get_alignment_path_prefix_by_ResultSet($rset);

  my $bam_file = $path . '.bam';
  my $bed_file = $path . '.bed';
  
  if(! -e $bam_file) {
    confess("Can't find bam file $bam_file!");
  }
  
  if (! -e $bed_file) {
    my $cmd = qq(bamToBed -i $bam_file > ${bed_file}.part);
    run_system_cmd($cmd);
    $cmd = qq(mv ${bed_file}.part $bed_file);
    run_system_cmd($cmd);
    push @file_to_delete_after_cell_line_has_been_processed, $bed_file;
  }
  
  push @file_to_delete_after_cell_line_has_been_processed, $bed_file;
  
  foreach my $current_file (@file_to_delete_after_cell_line_has_been_processed) {
    $self->dataflow_output_id( {
	file_to_delete_after_cell_line_has_been_processed => $current_file,
    }, 7);
  }
  
  my $align_files = $self->get_alignment_files_by_ResultSet_formats($rset, ['bam']);
  $align_files->{bed} = $bed_file;

  $self->set_param_method('bam_file', $align_files->{bam}, 'required');  # For bai test/creation 

  $self->new_Importer_params( 
   {#-prepared            => 1, #Will be if derived from BAM in get_alignment_file_by_InputSets
    #but we aren't extracting the slice names in this process yet
    -input_feature_class => 'result', #todo remove this requirement as we have output_set?
    -format              => 'SEQUENCING',#This needs changing to different types of seq
    -recover             => 1, 
    -force               => 1, #for store_window_bins_by_Slice_Parser
    -parser              => 'Bed',
    -output_set          => $rset,
    -input_files         => [$align_files->{bed}], #Can we set these after init?
    -slices              => generate_slices_from_names
                             ($self->out_db->dnadb->get_SliceAdaptor, 
                              $self->slices, 
                              $self->skip_slices, 
                              'toplevel', 0, 1),#nonref, incdups
   });
  
  #todo why do we have to set EFG_DATA?
  #$ENV{EFG_DATA} = $self->output_dir;
  
  #This picks up the defaults param from the input_id
  my $imp = $self->get_Importer;
  
  if((! $self->param_silent('merge')) &&
    ($imp->prepared)){
    $self->param_required('slices');     
  }
  
  return;
}


sub run {   # Check parameters and do appropriate database/file operations... 
  my $self = shift; 

  my $Imp  = $self->get_Importer;
  
  if($self->param_silent('merge')){
    $self->throw_no_retry('CollectionWriter no longer supports collection merge mode');
  }
  if($Imp->prepared ){
    $self->throw_no_retry('CollectionWriter no longer supports collection generation');
  }
  #Prepare or write slice col
  
  my $expected_bam_index_file = $self->bam_file.'.bai';
  
  # Test/create bai index file
  if(! -e $expected_bam_index_file) {
    my $cmd = 'samtools index '.$self->bam_file;  # -b option not require for this version
    run_system_cmd($cmd);
  }
  $self->dataflow_output_id( {
      file_to_delete_after_cell_line_has_been_processed => $expected_bam_index_file,
  }, 7);

  $Imp->read_and_import_data('prepare');

  # This now only rebokes the IMPORTED status
  # so we will need to manage that before we can remove the Importer usage
  # We're not actually storing anything yet, so that can be done in
  # the PeakCaller/BigWigWriter?
  # This was also creating the prepared bed for collection generation
  # Looks like this wasn't being used for CCAT peak calling
  # There small potential that an unsorted bed file could be used for CCAT
  # Although this is always generated from the sorted bam at present,
  # as opposed to a pre-computed unverified/sorted bed file.


  my $output_id = {%{$self->batch_params}, 
		    # These are already param_required by fetch_Set_input
		    dbID         => $self->param('dbID'),
		    set_name     => $self->param('set_name'),  # mainly for readability
		    set_type     => $self->param('set_type'),
		    filter_from_format => undef,                 
		  }; 

  #Need to add a final semaphored job, to do the clean up
  #bed files only, as we keep the bams
  $self->branch_job_group(2, [$output_id]); #BigWigWriter data flow

  my $fset = $self->FeatureSet;
  
  if(defined $fset){
    $self->branch_job_group('run_'.$fset->analysis->logic_name, [{%$output_id}]);
  }
  
  return;
}


sub write_output { 
  my $self = shift;    
  $self->dataflow_job_groups;
  return; 
}



1;
