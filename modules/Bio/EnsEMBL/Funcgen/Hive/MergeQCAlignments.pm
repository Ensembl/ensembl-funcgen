=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Hive::Funcgen::MergeQCAlignments

=head1 DESCRIPTION

Merges bam alignments from split (replicate or merged) fastq files.

=cut

package Bio::EnsEMBL::Funcgen::Hive::MergeQCAlignments;

use warnings;
use strict;

use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd );
use Bio::EnsEMBL::Funcgen::Sequencing::SeqTools;# merge_bams 

use base qw( Bio::EnsEMBL::Funcgen::Hive::BaseDB );

my %valid_flow_modes = (
  replicate => undef,
  merged    => undef,
  signal    => undef
); 

sub fetch_input {  
  my $self = shift;
  
  $self->SUPER::fetch_input();
  my $result_set = $self->fetch_Set_input('ResultSet');
  
  $self->get_param_method('bam_files',  'silent');
  $self->get_param_method('fastq_files',  'silent');
    
  if((! $self->bam_files ) && $self->fastq_files) {
    my $bam_file;
    $self->bam_files([ map {($bam_file = $_) =~ s/\.fastq_([0-9]+)$/.$1.bam/o; $bam_file} @{$self->fastq_files} ]);
  } elsif(! $self->bam_files) {
    $self->throw_no_retry('No bam_files or fastq_files have been defined');    
  }
  
  $self->get_param_method('set_prefix', 'required');  #This is control specific
  my $flow_mode = $self->get_param_method('flow_mode',  'required');
  $self->set_param_method('run_controls', 0); 
  
  if(! exists $valid_flow_modes{$flow_mode}) {
    throw("$flow_mode is now a valid flow_mode parameter, please specify one of:\t".
      join(' ', keys %valid_flow_modes));  
  }
  
  if(($flow_mode ne 'signal') && $self->get_param_method('result_set_groups', 'silent')) {
    throw("The $flow_mode flow mode is not valid for use with result_set_groups"); 
  }

  if($flow_mode eq 'replicate') {
    $self->get_param_method('permissive_peaks', 'required');
  }
 
  $self->init_branching_by_analysis;
  return;
}


sub run {
  my $self       = shift;
  my $result_set = $self->ResultSet;
  my $cmd;

  return;
  
  ### CLEAN FASTQS ###
  if($self->fastq_files){
    #Run with no exit flag so we don't fail on retry
    $cmd = 'rm -f '.join(' ', @{$self->fastq_files});
    $self->helper->debug(3, "Removing fastq chunks:\n$cmd");
    run_system_cmd($cmd, 1);
  }
  
  ### MERGE BAMS ###
  my $file_prefix  = $self->get_alignment_path_prefix_by_ResultSet($result_set, $self->run_controls); 
  my $unfiltered_bam     = $file_prefix.'.bam';
  
  $self->helper->debug(1, "Merging bams to:\t".$unfiltered_bam); 
  
  merge_bams({
    input_bams => $self->bam_files, 
    output_bam => $unfiltered_bam, 
    debug => $self->debug,
  });
  
  my $tmp_bam = "${unfiltered_bam}.nodups.bam";
  
  remove_duplicates_from_bam({
    input_bam  => $unfiltered_bam,
    output_bam => $tmp_bam, 
    debug      => $self->debug,
  });
  
  unlink($unfiltered_bam);
  run_system_cmd("mv $tmp_bam $unfiltered_bam");

  return;
}

sub write_output {

  my $self = shift;

  my $result_set   = $self->ResultSet;
  my %batch_params = %{$self->batch_params};
  my $flow_mode    = $self->flow_mode;
  
  if ($flow_mode eq 'merged') {
  
    my %output_id = (
      set_type    => 'ResultSet',
      set_name    => $self->ResultSet->name,
      dbID        => $self->ResultSet->dbID,
      garbage     => $self->bam_files,
    );

    $self->branch_job_group('DefineMergedDataSet', [{%batch_params, %output_id}]);
  }
  
  if ($flow_mode eq 'replicate') {
  
    my %output_id = (
      set_type      => 'ResultSet',
      set_name      => $self->ResultSet->name,
      dbID          => $self->ResultSet->dbID,
      garbage       => $self->bam_files,
      
      # If flowing to replicate processing, set the permissive_peaks parameter.
      # permissive_peaks is always set to "SWEmbl_R0005" in our runs.
      #
      peak_analysis => $self->permissive_peaks,
    );

    $self->branch_job_group('run_'.$self->permissive_peaks.'_replicate', [{%batch_params, %output_id}]);
  }
  
  if ($flow_mode eq 'signal') {
    # Should autoflow the input, creation of jobs is now done in 
    # Bio::EnsEMBL::Funcgen::Hive::JobFactorySignalProcessing
    #
    return;
  }
  
  $self->dataflow_job_groups;
  return;
}

1;
