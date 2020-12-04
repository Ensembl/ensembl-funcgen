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

package Bio::EnsEMBL::Funcgen::ExampleFeature;

use strict;
use Bio::EnsEMBL::Utils::Exception qw( deprecate );
use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_get_or_set
  _generic_fetch
);

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

sub _constructor_parameters {
  return {
    dbID                => 'dbID',
    adaptor             => 'adaptor',
    ensembl_object_type => 'ensembl_object_type',
    ensembl_id          => 'ensembl_id',
  };
}

sub dbID                { return shift->_generic_get_or_set('dbID', @_); }
sub adaptor {return shift->_generic_get_or_set('adaptor', @_);}
sub ensembl_object_type { return shift->_generic_get_or_set('ensembl_object_type', @_); }
sub ensembl_id          { return shift->_generic_get_or_set('ensembl_id', @_); }


#=head2 summary_as_hash
#
#  Example       : $summary = $frip->summary_as_hash;
#  Description   : Returns summary in a hash reference.
#  Returns       : Hashref of descriptive strings
#  Status        : Intended for internal use (REST)
#
#=cut
#
#sub summary_as_hash {
#  my $self   = shift;
#  
#  my $peak_calling = $self->fetch_PeakCalling;
#  
#  return {
#    ensembl_object_type  => $self->ensembl_object_type,
#    ensembl_ids         => $ensembl_ids->summary_as_hash,
#  };
#}


1;
