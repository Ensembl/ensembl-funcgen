
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

Bio::EnsEMBL::Funcgen::Hive::PostprocessIDR

=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Funcgen::Hive::PostprocessIDR;

use warnings;
use strict;
 
use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( scalars_to_objects );
use base ('Bio::EnsEMBL::Funcgen::Hive::BaseDB');

#TODO 
# 1 ADD IDR based QC
# 2 Build output id

sub fetch_input {   # fetch parameters...
  my $self = shift;
  $self->SUPER::fetch_input;
  
  #From SubmitIDR
  #{%{$batch_params},
  #                            dbIDs    => \@rset_ids,
  #                            set_type => 'ResultSet'}

  #Can either get the accu idr peaks or read from file
  #based on the 
    
  #How are we going to differentiate pseudo reps?
  #SubmitIDR will know about these,as they will be passed explicitly from GeneratePseudoReplicates
  #as file names (not sets!?).
  
  
  if($self->param_required('set_type') ne 'ResultSet'){
    throw("Param set_type should be 'ResultSet', not:\t".$self->param('set_type'));  
  }
  
  my $fset_ids = $self->get_param_method('dbIDs',  'required');
  assert_ref($fset_ids, 'ARRAY', 'FeatureSet dbIDs');
  $self->get_param_method('permissive_peaks', 'required');
  $self->get_param_method('idr_peaks',        'required');
  my $idr_peaks = $self->get_param_method('idr_peak_counts', 'required'); 
  #This is accumulated data from the RunIDR fan jobs, submitted & semaphored from PreprocessIDR
  assert_ref($idr_peaks, 'ARRAY', 'IDR peaks');
 
  return;
}


sub run {   # Check parameters and do appropriate database/file operations... 
  my $self  = shift;  
  my $batch_params  = $self->batch_params;
  my $rsets         = &scalars_to_objects($self->out_db,  'ResultSet',
                                          'fetch_by_dbID', $self->dbIDs);
  my $peak_analysis = &scalars_to_objects($self->out_db, 'Analysis',
                                          'fetch_by_logic_name',
                                          $self->permissive_peaks)->[0];                                                

  my $max_peaks = (sort {$a <=> $b} @{$self->idr_peak_counts})[-1];#Take the highest!


  #Where is merging going to happen?
  #We should do this in DefineResultSet
  #Pooled file should already exist from Pseudo-rep creation
  #but handle absence, just incase we want to run without this step

  
  #todo IDR based QC here! 
  
  
  #Need to flow peak_analysis here, such that we over-ride the non-IDR defaults
  #If we leave as is, with run_idr enable, PreprocessAlignments will
  #use the pre-IDR analysis.
  
  #Is there any other way to handle this?
  #the run_idr param is now a bit ambiguous
  #Should this be turned into a string value pre|post?
  #Probably easier just to over-ride peak analysis defaults with peak_analysis param
  #
  
  #Do the alignment merge here as we already have access to the rep ResultSets here?
  #or do it in 
  
  my @controls = grep {$_->is_control} @{$rsets->[0]->get_support};
  #my (@signals, @rep_bams);  
  my @signals  = grep {! $_->is_control} (map {$_->get_support} @$rsets);
  
  #my $filter_format = $self->param_silent('bam_filtered') ? undef : 'bam'; 

  my $set_prefix = $self->get_set_prefix_from_Set($rsets->[0]); 

  #We can't get the alignment file from get_alignment_files_by_ResultSet_formats or 
  #get_alignment_file_prefix_by_ResultSet as the ResultSet doesn't exist yet

  #This is why we were doing the merge in DefineResultSets


  #Let's do it there instead...GRR!

  
  #foreach my $rep_rset(@$rsets){
  #  
  #  push @signals, grep {! $_->is_control} @{$rep_rset->get_support};  
  #  
  #  #This should only return 1 bam file else throw
  #  push @rep_bams, 
  #    @{$self->get_alignment_files_by_ResultSet_formats($rep_rset, ['bam'],
  #                                                      undef,  # control flag
  #                                                      undef,  # all_formats flag
  #                                                      $filter_format)};
  #}
  
  #This is flowing to DefineMergedReplicateResultSet
  $self->branch_job_group(2, {%{$batch_params},
                              max_peaks                  => $max_peaks, #This is batch flown
                              peak_analysis              => $self->idr_peaks,
                              merge_replicate_alignments => 1,
                             
                              input_subset_ids => 
                               {
                                controls     => \@controls,
                                $set_prefix  => \@signals, }});
                                 
  return;
}



sub write_output {  # Create the relevant jobs
  shift->dataflow_job_groups;
  return;
}

1;