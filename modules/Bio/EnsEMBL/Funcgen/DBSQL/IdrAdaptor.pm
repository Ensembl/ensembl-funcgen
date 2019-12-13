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

package Bio::EnsEMBL::Funcgen::DBSQL::IdrAdaptor;

use strict;
use base 'Bio::EnsEMBL::Funcgen::DBSQL::GenericAdaptor';

sub object_class {
    return 'Bio::EnsEMBL::Funcgen::Idr';
}

sub _tables {
  return ['idr', 'i']
}

sub fetch_all_failed {
  my $self         = shift;
  
  my $features = $self->fetch_all(
    "failed_idr_pairs is not null and max_peaks is null"
  );
  return $features;
}

sub _fetch_by_experiment_id {

  my $self          = shift;
  my $experiment_id = shift;
  
  my $feature = $self->fetch_all(
    "experiment_id = " . $experiment_id
  );
  return $feature->[0];
}

sub _fetch_by_Experiment {

  my $self       = shift;
  my $experiment = shift;
  
  return $self->_fetch_by_experiment_id($experiment->dbID);
}

1;
