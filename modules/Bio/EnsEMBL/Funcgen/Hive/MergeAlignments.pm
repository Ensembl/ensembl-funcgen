=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Hive::Funcgen::MergeAlignments

=head1 DESCRIPTION

Merges bam alignments from split (replicate or merged) fastq files.

=cut

package Bio::EnsEMBL::Funcgen::Hive::MergeAlignments;

use strict;

use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd );
use Bio::EnsEMBL::Funcgen::Sequencing::SeqTools;# merge_bams 

use base qw( Bio::EnsEMBL::Funcgen::Hive::BaseDB );

sub fetch_input {

  my $self = shift;

  $self->SUPER::fetch_input();
  my $result_set = $self->fetch_Set_input('ResultSet');

  $self->get_param_method('bam_files',  'silent');
  $self->get_param_method('fastq_files',  'silent');

  if((! $self->bam_files ) && $self->fastq_files) {
    my $bam_file;
    
    # Creates the name of the bam files from the names of the fastq files
    #
    $self->bam_files([ map {($bam_file = $_) =~ s/\.fastq_([0-9]+)$/.$1.bam/o; $bam_file} @{$self->fastq_files} ]);
    
  } elsif(! $self->bam_files) {
    $self->throw_no_retry('No bam_files or fastq_files have been defined');
  }
  $self->get_param_method('run_controls',  'required');
  return;
}


sub run {
  my $self       = shift;
  my $result_set = $self->ResultSet;
  my $cmd;
  
  if($self->fastq_files) {
  
    # Run with no exit flag so we don't fail on retry
    $cmd = 'rm -f '.join(' ', @{$self->fastq_files});
    $self->helper->debug(3, "Removing fastq chunks:\n$cmd");
    run_system_cmd($cmd, 1);
  }
  
  my $bam_file_with_unmapped_reads_and_duplicates = $self->param('bam_file_with_unmapped_reads_and_duplicates');
  
  $self->helper->debug(1, "Merging bams to:\t".$bam_file_with_unmapped_reads_and_duplicates); 

  # HACK: Only merge the files, if the merged file doesn't already exists.
  #
  # This can happen, because experiments can share a control. If one of these
  # experiments was run, then the bam file for the control has already been
  # generated.
  #
  if (! -e $bam_file_with_unmapped_reads_and_duplicates) {
    merge_bams({
      input_bams => $self->bam_files, 
      output_bam => $bam_file_with_unmapped_reads_and_duplicates,
      debug => $self->debug,
    });
  }
  

#   my $cmd=qq(java picard.cmdline.PicardCommandLine ValidateSamFile INPUT=$bam_file_with_unmapped_reads_and_duplicates);
#   $self->hive_run_system_cmd($cmd);
  
  foreach my $current_bam_file (@{$self->bam_files}) {
    unlink($current_bam_file);
  }

  return;
}

1;
