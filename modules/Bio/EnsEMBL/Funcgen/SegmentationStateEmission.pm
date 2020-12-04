=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

=head1 SYNOPSIS

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::SegmentationStateEmission;

use strict;
use Bio::EnsEMBL::Utils::Exception qw( deprecate );
use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_get_or_set
);

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

sub _constructor_parameters {
  return {
    dbID     => 'dbID',
    adaptor  => 'adaptor',
    state    => 'state',
    segmentation => 'segmentation',
    CTCF     => 'CTCF',
    DNase1   => 'DNase1',
    H3K27ac  => 'H3K27ac',
    H3K27me3 => 'H3K27me3',
    H3K36me3 => 'H3K36me3',
    H3K4me1  => 'H3K4me1',
    H3K4me2  => 'H3K4me2',
    H3K4me3  => 'H3K4me3',
    H3K9ac   => 'H3K9ac',
    H3K9me3  => 'H3K9me3',
  };
}

sub dbID     { return shift->_generic_get_or_set('dbID',     @_); }
sub adaptor {return shift->_generic_get_or_set('adaptor', @_);}
sub state        { return shift->_generic_get_or_set('state',    @_); }
sub segmentation { return shift->_generic_get_or_set('segmentation',    @_); }
sub CTCF     { return shift->_generic_get_or_set('CTCF',     @_); }
sub DNase1   { return shift->_generic_get_or_set('DNase1',   @_); }
sub H3K27ac  { return shift->_generic_get_or_set('H3K27ac',  @_); }
sub H3K27me3 { return shift->_generic_get_or_set('H3K27me3', @_); }
sub H3K36me3 { return shift->_generic_get_or_set('H3K36me3', @_); }
sub H3K4me1  { return shift->_generic_get_or_set('H3K4me1',  @_); }
sub H3K4me2  { return shift->_generic_get_or_set('H3K4me2',  @_); }
sub H3K4me3  { return shift->_generic_get_or_set('H3K4me3',  @_); }
sub H3K9ac   { return shift->_generic_get_or_set('H3K9ac',   @_); }
sub H3K9me3  { return shift->_generic_get_or_set('H3K9me3',  @_); }

sub get_SegmentationStateAssignment {
    my $self = shift;

    my $segmentation_state_assignment_adaptor =
        $self->adaptor->db->get_SegmentationStateAssignmentAdaptor;
    if (! defined $segmentation_state_assignment_adaptor) {
        throw("Couldn't get an SegmentationStateAssignmentAdaptor!");
    }
    my $segmentation_state_assignment = $segmentation_state_assignment_adaptor
    ->fetch_by_state_segmentation(
        $self->state,
        $self->segmentation,
    );

    return $segmentation_state_assignment;
}

1;
