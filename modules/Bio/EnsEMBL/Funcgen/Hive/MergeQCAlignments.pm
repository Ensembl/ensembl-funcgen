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
  <http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Hive::Funcgen::MergeQCAlignments

=head1 DESCRIPTION



=cut

package Bio::EnsEMBL::Funcgen::Hive::MergeQCAlignments;

use warnings;
use strict;

use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd );

use base qw( Bio::EnsEMBL::Funcgen::Hive::BaseDB );

#TODO... use and update the tracking database dependant on no_tracking...

#todo
# Reimplement repository support
# Status handling/setting?


sub fetch_input {  
  my $self = shift;
  $self->SUPER::fetch_input();
  my $rset = $self->fetch_Set_input('ResultSet');
  $self->get_param_method('output_dir', 'required');
  $self->get_param_method('bam_files',  'required');
  $self->get_param_method('set_prefix', 'required');#was chunk_fiule_prefix in MergeChunkResultSetFastqs
  my $flow_mode = $self->get_param_method('flow_mode',  'required');
  #could have recreated output_dir and merged_file_name from ResultSet and run_controls
  #as we did in MergeChunkResultSetFastqs, but passed for convenience

  my $species = $self->species; #set in _set_outdb
  my $gender  = $rset->cell_type->gender || 'male';
  
  my $sam_header = join('/', ($self->param_required('data_root_dir'),
                               'sam_header',
                               $species,
                               $species.'_'.$gender.'_'.$self->assembly.'_unmasked.header.sam'));
  $self->set_param_method('sam_header', $sam_header);

  if($flow_mode ne 'signal'){
    $self->get_param_method('result_set_groups', 'required');
  }
  
 
  #my $repository = $self->_repository();
  #if(! -d $repository){
  #  system("mkdir -p $repository") && throw("Couldn't create directory $repository");
  #}

    
  return;
}

#Is the cating of rep numbers going to cause problems?
#The Peaks and Collections pipelines need to be able to 
#recreate these file names.
#in fact we need to share the same code used by get_alignment_file_by_ResultSet_formats
#which does not use the rep numbers
#This may be risky
#In the event of a disconnect between update of an ResultSet with new reps
#and running the peaks/collections on a pre-exitingalignment file(without new rep)
#This could also be managed with a ResultSet status. Any updates to the ResultSet
#would require unsetting the ALIGNED status

sub run {
  my $self         = shift;
  my $rset         = $self->ResultSet;
  my $sam_header   = $self->sam_header;
  my @bam_files    = @{$self->bam_files};  
  my $bam_file     = $self->output_dir.'/'.$self->set_prefix.'.unfiltered.bam';
  my $filtered_bam = $self->output_dir.'/'.$self->set_prefix.'.bam';
  run_system_cmd("samtools merge -f -h $sam_header $bam_file @bam_files");
  #Maybe remove duplicates as default behavior ?.
  #my $merge_cmd="samtools merge -h $sam_header - ${file_prefix}.[0-9]*.sorted.bam | ";
  #$merge_cmd.=" samtools rmdup -s - ${file_prefix}.sorted.bam";


  #keep end result as bam... it's more compact: only when passing to the peak caller we pass to sam or other format...
  #merge_cmd.=" samtools view -h - | gzip -c > ${file_prefix}${align_type}.sam.gz"


#todo convert this to wite to a result_set_report table

  my $alignment_log = $self->output_dir.'/'.$self->set_prefix.".alignment.log";

  my $cmd="echo \"Alignment QC - total reads as input: \" >> ${alignment_log}".
    ";samtools flagstat $bam_file | head -n 1 >> ${alignment_log}".
    ";echo \"Alignment QC - mapped reads: \" >> ${alignment_log} ".
    ";samtools view -u -F 4 $bam_file | samtools flagstat - | head -n 1 >> ${alignment_log}".
    "; echo \"Alignment QC - reliably aligned reads (mapping quality >= 1): \" >> ${alignment_log}".
    ";samtools view -u -F 4 -q 1 $bam_file | samtools flagstat - | head -n 1 >> ${alignment_log}";
  #Maybe do some percentages?
  run_system_cmd($cmd);
  
  #my $repository = $self->_repository();
  #move("${file_prefix}.sorted.bam","${repository}/${set_name}.samse.bam");
  #my $convert_cmd =  "samtools view -h ${file_prefix}.sorted.bam | gzip -c - > ${repository}/${set_name}.samse.sam.gz";

  run_system_cmd("rm -f @bam_files");


  #Filter and QC here FastQC or FASTX?
  #filter for MAPQ >= 15 here? before or after QC?
  #PhantomPeakQualityTools? Use estimate of fragment length in the peak calling?

  warn "Need to implement post alignment QC here. Filter out MAPQ <16. FastQC/FASTX/IDR?/PhantomPeakQualityTools frag length?";
  #Add ResultSet status setting here!
  #ALIGNED
  #ALIGNMENT_QC_OKAY





  my $flow_mode    = $self->flow_mode;
  my %batch_params = %{$self->batch_params};
  
  if($flow_mode ne 'signal'){
    #Data flow to DefineDataSet.pm i.e. DefineMergedOutputSet or DefineReplicateOutputSet
    $self->branch_job_group('Define'.$flow_mode.'DataSet', 
                            [{%batch_params,
                              set_type    => 'ResultSet',
                              set_name    => $self->ResultSet->name,
                              dbID        => $self->ResultSet->dbID}]);
  }
  else{ #Run signal fastqs
    my $rset_groups = $self->result_set_groups;
    my $align_anal  = $rset->analysis->logic_name;
    
    #This bit is very similar to some of the code in DefineResultSets
    #jsut different enough not to be subbable tho 

    foreach my $rset_group(keys %{$rset_groups}){
      my @fan_jobs;
        
      for my $i(0...$#{$rset_groups->{$rset_group}{dbIDs}}){
        push @fan_jobs, {%batch_params,
                         set_type    => 'ResultSet',
                         set_name    => $rset_groups->{$rset_group}{set_names}->[$i],
                         dbID        => $rset_groups->{$rset_group}{dbIDs}->[$i]};
      }
         
      my $branch = ($rset_group eq 'merged') ? 
        'Preprocess_'.$align_anal.'_merged' : 'Preprocess_'.$align_anal.'_replicate';   
      
      my @job_group = ($branch, \@fan_jobs);   
         
      if($rset_group ne 'merged'){ #Add RunIDR semaphore
        push @job_group, ('RunIDR', 
                          [{%batch_params,
                            dbIDs     => $rset_groups->{$rset_group}{dbIDs},
                            set_names => $rset_groups->{$rset_group}{set_names},
                            set_type  => 'ResultSet'}]);     
      }  
       
      $self->branch_job_group(@job_group);
    }
  }

  return;
}


sub write_output {  # Create the relevant jobs
  shift->dataflow_job_groups;
  return;
}



1;
