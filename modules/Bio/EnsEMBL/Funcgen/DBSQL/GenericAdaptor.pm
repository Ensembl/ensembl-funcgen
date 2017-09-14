=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Funcgen::DBSQL::GenericAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw );
use base 'Bio::EnsEMBL::Funcgen::DBSQL::BaseAdaptor';

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::DBSQL::GenericAdaptorMethods';

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  $self->init_generic_adaptor(@args);
  return $self;
}

sub object_class {
  my $self = shift;
  throw('object_class has to be overwritten by the subclass in ('. ref $self .')!');
}

1;
