
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
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( scalars_to_objects );
use base ('Bio::EnsEMBL::Funcgen::Hive::BaseDB');


sub fetch_input {   # fetch parameters...
  my $self = shift;
  $self->SUPER::fetch_input;
  
  
  #Hmm, peaks should already have been called for each of these result sets
  #so we can simply bring back the data/feature sets, and run IDR on them
  #{%batch_params,
  #                          dbIDs     => $rset_groups->{$rset_group}{dbIDs},
  #                          set_names => $rset_groups->{$rset_group}{set_names},
  #                          set_type  => 'ResultSet'}
  if($self->param_required('set_type') ne 'ResultSet'){
    throw('');  
  }
  
  my $rset_ids = $self->get_param_method('dbIDs',  'required');
  assert_ref($rset_ids, 'ARRAY', 'ResultSet dbIDs');
  
  #IDR analysis should probably be specified as a default/pipeline_wide and batch flowable logic_name
  #This will need to be defined in BaseSequenceANalysis as it is needed by several confs
  
  $self->set_param_method('idr_analysis', 'required');
      


  return;
}


sub run {   # Check parameters and do appropriate database/file operations... 
  my $self         = shift;
  #my $helper    = $self->helper;
  my $idr_analysis = &scalars_to_objects($self->out_db, 'Analysis',
                                                        'fetch_by_logic_name',
                                                        [$self->idr_analysis])->[0];
  
  my $rsets = &scalars_to_objects($self->out_db, 'ResultSet',
                                                 'fetch_by_dbID',
                                                 $self->dbIDs);
    if(! &_are_controls($ctrls)){
      throw("Found unexpected non-control InputSubsets specified as controls\n\t".
        join("\n\t", map($_->name, @$ctrls)));
    }
  }
                                             
  #Get the FeatureSets for each ResultSet.
  my @fsets;
  my $dset_a = $self->out_db->get_DataSetAdaptor;
  my $throw = '';
  
  foreach my $rset(@rsets){
    my @dsets = @{$dset_a->fetch_all_by_ResultSet};
    
    if( scalar(@dsets) != 1 ){
      $throw .= "Could not find unqiue DataSet assoicated to ResultSet:\t".$rset->name."\n";
    }
    
    my $fset = $dsets[0]->product_FeatureSet;
    
    if(! defined $fset){
      $throw .= "Could not find associated FeatureSet for ResultSet:\t".$rset->name."\n";  
    }
    else{
      push @fsets, $fset;  
    }
  }
    
  if($throw){
    throw($throw);  
  }  
    
  
  #Now prepare the input for the IDR analysis
  #Dump to bed?
  
  
    #        $self->branch_job_group($branch, [{%{$batch_params},
    #                                           dbID       => [$rset->dbID], 
    #                                           set_name   => [$rset->name],
    #                                           set_type    => 'ResultSet'}]);){
      
  return;
}



sub write_output {  # Create the relevant jobs
  shift->dataflow_job_groups;
  return;
}

1;