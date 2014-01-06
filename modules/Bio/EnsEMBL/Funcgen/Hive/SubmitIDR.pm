
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

This is a very simple module which handles IDR job submissions and enables the 
correct semaphore config (fan/funnel). 

# This could be omitted in multi-sempahores are supported e.g.
#    '2->A'    => ['DefineReplicateDataSet'], #fan
#    'A->3->B' => ['RunIDR'],                 #fan and immediate funnel
#    'B->4'    => ['PostProcessIDRReplicates],

=cut

package Bio::EnsEMBL::Funcgen::Hive::SubmitIDR;

use warnings;
use strict;
 
use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( scalars_to_objects );
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
                                                               $self->param_required('permissive_peaks'));
  return;
}


sub run {   # Check parameters and do appropriate database/file operations... 
  my $self  = shift;  
  my $rsets = &scalars_to_objects($self->out_db, 'ResultSet',
                                                 'fetch_by_dbID',
                                                 $self->dbIDs);
  my $peak_analysis = $self->peak_analysis;                                                 
  my $dset_adaptor = $self->out_db->get_DataSetAdaptor;
  my $batch_params = $self->batch_params;
  my @fset_ids;

  foreach my $rset(@rsets){
    
    my $dsets = $dset_adaptor->fetch_all_by_supporting_set($rset);
    #We need to peak analysis here to be able to support >1 data/feature set
    #Currently this is a pipelinewide permissive_peaks param
    #so would need to make this batch flowable, to allow different permissive analyses       
    my $found_fset = 0;
    
    if(scalar(@$dsets) == 0){
      throw("Could not find an associated DataSet for ResultSet:\t".$rset->name);  
    }
    
    
    foreach my $dset(@$dsets){
      my $fset = $dset->get_product_FeatureSet;
      
      if($fset->analysis->dbID == $peak_analysis->dbID){
        $found_fset = 1;
        last;
      }
    }

    if(! $found_fset){
      throw('Could not find '.$peak_analysis->logic_name." Data/FeatureSet for ResultSet:\t".$rset->name);  
    }    
  
    
    #TODO check the feature_set_stat or statuses
    #such that we know the peak calling jobs has finished and passed QC!    
    #Do this for each before we submit IDR jobs, as we may have to drop some reps
    push @fset_ids, $fset->dbID;
  }


  #Build 2 way rep combinations for IDR jobs
  my @idr_job_ids;
  my $last_i = $#fset_ids - 1;
  
  foreach my $i(0..$last_i){
    my $next_i = $i + 1;
  
    foreach my $j($next_i..$#fset_ids){
      push @idr_job_ids, {{%{$batch_params},
                          dbIDs    => [$fset_ids->[i], $fset_ids->[$j]],
                          set_type => 'FeatureSet'};
    }  
  }
  
  
  $self->branch_job_group(2, \@idr_job_ids,
                          3, {%{$batch_params},
                              dbIDs    => \@fset_ids,
                              set_type => 'FeatureType'});
      
  return;
}



sub write_output {  # Create the relevant jobs
  shift->dataflow_job_groups;
  return;
}

1;