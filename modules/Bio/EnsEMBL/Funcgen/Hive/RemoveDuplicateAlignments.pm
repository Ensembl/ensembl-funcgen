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

Bio::EnsEMBL::Hive::Funcgen::RemoveDuplicateAlignments

=cut

package Bio::EnsEMBL::Funcgen::Hive::RemoveDuplicateAlignments;

use strict;

use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd );
use Bio::EnsEMBL::Funcgen::Sequencing::SeqTools;

use base qw( Bio::EnsEMBL::Funcgen::Hive::BaseDB );

sub fetch_input {

  my $self = shift;
  $self->SUPER::fetch_input();  
  my $result_set = $self->fetch_Set_input('ResultSet');
  $self->get_param_method('run_controls',  'required');
  return;
}


sub run {
  my $self       = shift;
  my $result_set = $self->ResultSet;
  my $cmd;
  
  my $file_prefix  = $self->get_alignment_path_prefix_by_ResultSet($result_set, $self->run_controls); 
  my $unfiltered_bam     = $file_prefix.'.bam';
  
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

1;
