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

package Bio::EnsEMBL::Funcgen::PeakCallingStatistic;

use strict;
use Bio::EnsEMBL::Utils::Exception qw( deprecate );
use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_get_or_set
);

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

sub _constructor_parameters {
  return {
    dbID            => 'dbID',
    adaptor         => 'adaptor',
    peak_calling_id => 'peak_calling_id',
    epigenome_id    => 'epigenome_id',
    feature_type_id => 'feature_type_id',
    statistic       => 'statistic',
    value           => 'value',
  };
}

sub dbID            { return shift->_generic_get_or_set('dbID',            @_); }
sub adaptor {return shift->_generic_get_or_set('adaptor', @_);}
sub peak_calling_id { return shift->_generic_get_or_set('peak_calling_id', @_); }
sub statistic       { return shift->_generic_get_or_set('statistic',       @_); }
sub value           { return shift->_generic_get_or_set('value',           @_); }
sub epigenome_id    { return shift->_generic_get_or_set('epigenome_id',    @_); }
sub feature_type_id { return shift->_generic_get_or_set('feature_type_id', @_); }

sub get_PeakCalling {
  my $self = shift;
  
  my $peak_calling_adaptor = $self->adaptor->db->get_PeakCallingAdaptor;
  if (! defined $peak_calling_adaptor) {
    throw("Couldn't get an PeakCallingAdaptor!");
  }
  my $peak_calling = $peak_calling_adaptor->fetch_by_PeakCallingStatistic($self);
  return $peak_calling;
}

1;
