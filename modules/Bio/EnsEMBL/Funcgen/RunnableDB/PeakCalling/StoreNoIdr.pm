
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

Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::StoreNoIdr

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::StoreNoIdr;

use warnings;
use strict;
use base 'Bio::EnsEMBL::Hive::Process';

use Data::Dumper;

sub run {
  my $self = shift;
  my $species        = $self->param_required('species');
  my $execution_plan = $self->param_required('execution_plan');
  
  my $experiment_name = $execution_plan->{idr}->{name};
  
  my $experiment_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'experiment');
  my $experiment = $experiment_adaptor->fetch_by_name($experiment_name);
  if (! defined $experiment) {
    die;
  }
  
  use Bio::EnsEMBL::Funcgen::Idr qw ( NO_IDR  );
  my $idr = Bio::EnsEMBL::Funcgen::Idr->new(
    -experiment_id => $experiment->dbID,
    -max_peaks     => undef,
    -type          => NO_IDR,
  );
  my $idr_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'Idr');
  $idr_adaptor->store($idr);
  return;
}

1;
