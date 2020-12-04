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

package Bio::EnsEMBL::Funcgen::ProbeMapping;

use strict;
use Bio::EnsEMBL::Utils::Exception qw( deprecate );
use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_get_or_set
);

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

sub _constructor_parameters {
  return {
    dbID                => 'dbID',
    adaptor             => 'adaptor',
    gene_build_version  => 'gene_build_version',
    five_prime_utr      => 'five_prime_utr',
    three_prime_utr     => 'three_prime_utr',
    sample_probe_id     => 'sample_probe_id',
    sample_probe_set_id => 'sample_probe_set_id',
    release_version     => 'release_version',
    release_date        => 'release_date',
    assembly            => 'assembly',
  };
}

sub dbID                { return shift->_generic_get_or_set('dbID',                @_); }
sub adaptor {return shift->_generic_get_or_set('adaptor', @_);}
sub gene_build_version  { return shift->_generic_get_or_set('gene_build_version',  @_); }
sub five_prime_utr      { return shift->_generic_get_or_set('five_prime_utr',      @_); }
sub three_prime_utr     { return shift->_generic_get_or_set('three_prime_utr',     @_); }
sub sample_probe_id     { return shift->_generic_get_or_set('sample_probe_id',     @_); }
sub sample_probe_set_id { return shift->_generic_get_or_set('sample_probe_set_id', @_); }
sub release_version     { return shift->_generic_get_or_set('release_version',     @_); }
sub release_date        { return shift->_generic_get_or_set('release_date',        @_); }
sub assembly            { return shift->_generic_get_or_set('assembly',            @_); }

1;
