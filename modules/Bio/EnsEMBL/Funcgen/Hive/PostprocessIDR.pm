
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
  assert_ref($rset_ids, 'ARRAY', 'FeatureSet dbIDs');
  $self->get_param_method('permissive_peaks', 'required');
  my $idr_peaks = $self->get_param_method('idr_peaks', 'required'); 
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
                                          $self->param_required('permissive_peaks'))->[0];                                                

  my $max_peaks = (sort {$a <=> $b} @{$self->idr_peaks})[-1];#Take the highest!

  #Now build input job id for  DefineResultSet
  #Need to dataflow max_peaks
  #This needs to be batch flown from here
  #so we don't have to explicitly flow it in 
  #DefineMergedReplicateResultSet, DefineMergedDataSet and PreprocessAlignments

  #Where is merging going to happen?
  #We should do this in DefineResultSet
  #Pooled file should already exist from Pseudo-rep creation
  #but handle absence, just incase we want to run without this step

  
  #todo IDR based QC here! 
  
  
  
  $self->branch_job_group(2, {%{$batch_params},
                              max_peaks => $max_peaks});
      
  return;
}



sub write_output {  # Create the relevant jobs
  shift->dataflow_job_groups;
  return;
}

1;