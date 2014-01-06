
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
                                                 
#echo Format conversion 
# Note: you may wish to perlify these awk commands, beware that AWK numbering is 1-based
#awk 'BEGIN {OFS="\t"} NR == 9 {min=$5} NR > 14 {print $1,$2,$3,".",$7,".",($4-$8)/min,-1,-1,int($9-$1)}' $res1 > $res1_bed
#awk 'BEGIN {OFS="\t"} NR == 9 {min=$5} NR > 14 {print $1,$2,$3,".",$7,".",($4-$8)/min,-1,-1,int($9-$1)}' $res2 > $res2_bed

#This is picking up the min var from line 5, then processing from line 15 onwards.
#What is default SWEmbl output here?


#$res1 and $res2 here are the SWEmbl output. Do we depend on flat files here, or dump form the DB?
#

#    echo IDR analysis
#Rscript ~dz1/utils/idrCode/batch-consistency-analysis.r $res1_bed $res2_bed -1 $exp_name 0 F signal.value
#signal.value is ranking measure here                                                
    
    
#echo Compute cutoff value for IDR = 1%
# Same comment as above
#factor=`awk 'NF > 10 && $11 < 0.01 && (min == 0 || min > $5) {min = $5} END {print min}' ${exp_name}-overlapped-peaks.txt`
# This is when I realize that starting off in an actual scripting language would have been better... oh well
#new_cutoff="$(echo "$LOW_CUTOFF*$factor" | bc)"

#WTF Why is this subsampling? 
# echo Subsample
#samtools view -s 0.5 -b $bam1 > $sub_bam1
#samtools view -s 0.5 -b $bam2 > $sub_bam2
#samtools merge -f $merged_bam $sub_bam1 $sub_bam2
#rm -f $sub_bam1 $sub_bam2

#echo Final run with new -R param set to $new_cutoff
#SWEMBL -F -V -i $merged_bam -f 150 -R $new_cutoff -o $output -r $control 
 
    
    
#The following need to be done in the PostProcessIDRReplicates                                                 
#some batch consistency plots  
#Rscript batch-consistency-plot.r 3 /consistency/reps/chipSampleAllReps /consistency/reps/chipSampleRep1_VS_chipSampleRep2 /consistency/reps/chipSampleRep1_VS_chipSampleRep3 /consistency/reps/chipSampleRep2_VS_chipSampleRep3  


#
  
                                                 
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