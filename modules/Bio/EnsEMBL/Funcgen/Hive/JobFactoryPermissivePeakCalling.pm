=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

=cut

package Bio::EnsEMBL::Funcgen::Hive::JobFactoryPermissivePeakCalling;

use strict;

use base qw( Bio::EnsEMBL::Funcgen::Hive::BaseDB );

sub fetch_input {
  my $self = shift;
  
  # This sets out_db which is needed to get a ResultSetAdapter in the next 
  # command.
  #
  $self->SUPER::fetch_input();
  $self->fetch_Set_input('ResultSet');
  return;
}

sub run {
  my $self = shift;

  my $result_set   = $self->ResultSet;
  my %batch_params = %{$self->batch_params};

  my %output_id = (
    set_type      => 'ResultSet',
    set_name      => $self->ResultSet->name,
    dbID          => $self->ResultSet->dbID,
  );

  $self->branch_job_group(100, [{%batch_params, %output_id}]);
}

sub write_output {
  my $self = shift;
  $self->dataflow_job_groups;
  return;
}

1;
