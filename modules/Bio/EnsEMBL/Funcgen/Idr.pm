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

package Bio::EnsEMBL::Funcgen::Idr;

use strict;
use base qw( Exporter );
use vars qw( @EXPORT_OK );

our @EXPORT_OK = qw(

  IDR_ON_BIOLOGICAL_REPLICATES 
  IDR_ON_TECHNICAL_REPLICATES
  NO_IDR

);

use constant {

  IDR_ON_BIOLOGICAL_REPLICATES => 'on biological replicates',
  IDR_ON_TECHNICAL_REPLICATES  => 'on technical replicates',
  NO_IDR                       => 'no_idr',

  };
use Bio::EnsEMBL::Utils::Exception qw( deprecate );
use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_get_or_set
  _generic_fetch
);

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

sub _constructor_parameters {
  return {
    dbID           => 'dbID',
    adaptor        => 'adaptor',
    experiment_id => 'experiment_id',
    max_peaks     => 'max_peaks',
    type          => 'type',
    failed_idr_pairs => 'failed_idr_pairs',
  };
}

sub dbID          { return shift->_generic_get_or_set('dbID',          @_); }
sub adaptor {return shift->_generic_get_or_set('adaptor', @_);}
sub experiment_id { return shift->_generic_get_or_set('experiment_id', @_); }
sub max_peaks     { return shift->_generic_get_or_set('max_peaks',     @_); }
sub type          { return shift->_generic_get_or_set('type',          @_); }
sub failed_idr_pairs { return shift->_generic_get_or_set('failed_idr_pairs', @_); }

sub get_Experiment {
  return shift->_generic_fetch('experiment', 'get_ExperimentAdaptor', 'experiment_id');
}

=head2 summary_as_hash

  Example       : $summary = $idr->summary_as_hash;
  Description   : Returns summary in a hash reference.
  Returns       : Hashref of descriptive strings
  Status        : Intended for internal use (REST)

=cut

sub summary_as_hash {
  my $self   = shift;
  
  my $summary = {
    type      => $self->type,
  };
  
  if ($self->failed_idr_pairs) {
    $summary->{failed_idr_pairs} = $self->failed_idr_pairs;
  }
  if ($self->max_peaks) {
    $summary->{max_peaks} = $self->max_peaks;
  }
  return $summary;
}

1;
