=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Funcgen::DBSQL::ExampleFeatureAdaptor;

use strict;
use base 'Bio::EnsEMBL::Funcgen::DBSQL::GenericAdaptor';

sub object_class {
    return 'Bio::EnsEMBL::Funcgen::ExampleFeature';
}

sub _tables {
  return ['example_feature', 'ef']
}

sub fetch_all_by_ensembl_object_type {
    my ($self, $ensembl_object_type) = @_;
    
    my $constraint = 'ensembl_object_type  = ' . $ensembl_object_type;
    return $self->fetch_all($constraint);
}

1;
