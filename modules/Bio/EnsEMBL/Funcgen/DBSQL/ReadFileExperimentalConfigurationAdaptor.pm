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

Bio::EnsEMBL::DBSQL::Funcgen::RegulatoryBuildAdaptor

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::ReadFileExperimentalConfigurationAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use DBI qw(:sql_types);

use base 'Bio::EnsEMBL::DBSQL::BaseAdaptor';

sub _tables {
  return (
    ['read_file_experimental_configuration', 'rfec'  ],
  );
}

sub _columns {
  my $self = shift;
  
  return qw(
      rfec.read_file_experimental_configuration_id
      rfec.read_file_id
      rfec.experiment_id
      rfec.technical_replicate
      rfec.biological_replicate
  );
}

sub _default_where_clause {
  return '';
}

sub fetch_all_by_experiment_id {
  my $self = shift;
  my $experiment_id = shift;

  my $constraint = "rfec.experiment_id = ?";
  $self->bind_param_generic_fetch($experiment_id, SQL_INTEGER);
  
  my $object_list = $self->generic_fetch($constraint);
  
  if (!$object_list || @$object_list==0) {
    return;
  }
  return $object_list;
}

sub fetch_all_by_Experiment {
  my $self = shift;
  my $experiment = shift;
  my $experiment_id = $experiment->dbID;
  return $self->fetch_all_by_experiment_id($experiment_id);
}

sub _objs_from_sth {
  my ($self, $sth) = @_;

  my(
    $sth_fetched_dbID,
    $sth_fetched_read_file_id,
    $sth_fetched_technical_replicate,
    $sth_fetched_biological_replicate,
    $sth_fetched_experiment_id,
  );
  
  $sth->bind_columns (
    \$sth_fetched_dbID,
    \$sth_fetched_read_file_id,
    \$sth_fetched_experiment_id,
    \$sth_fetched_technical_replicate,
    \$sth_fetched_biological_replicate,
  );
  
  use Bio::EnsEMBL::Funcgen::ReadFileExperimentalConfiguration;
  
  my $experiment_adaptor = $self->db->get_ExperimentAdaptor();
  my $read_file_adaptor  = $self->db->get_ReadFileAdaptor();
  
  my @return_object_list;
  
  ROW: while ( $sth->fetch() ) {
  
    my $experiment = $experiment_adaptor->fetch_by_dbID($sth_fetched_experiment_id);
    my $read_file  = $read_file_adaptor->fetch_by_dbID($sth_fetched_read_file_id);

    my $current_object = Bio::EnsEMBL::Funcgen::ReadFileExperimentalConfiguration->new(
      -db                   => $self->db,
      -dbID                 => $sth_fetched_dbID,
      -read_file            => $read_file,
      -technical_replicate  => $sth_fetched_technical_replicate,
      -biological_replicate => $sth_fetched_biological_replicate,
      -experiment           => $experiment,
    );
    push @return_object_list, $current_object;
  }
  return \@return_object_list;
}

sub store {
  my ($self, @object) = @_;
  
  my $sth_store_object = $self->prepare("
    INSERT INTO read_file_experimental_configuration (
      read_file_id,
      experiment_id,
      biological_replicate,
      technical_replicate
    ) VALUES (?, ?, ?, ?)"
  );
  
  my $read_file_adaptor = $self->db->get_ReadFileAdaptor;
  
  foreach my $current_object (@object) {

    my $read_file    = $current_object->get_ReadFile;
    my $read_file_id = $read_file->dbID;
    
    if (! defined $read_file_id) {
    
      # If is has no id, then it hasn't been stored yet.
      #
      $read_file_adaptor->store($read_file);
      
      # Now it should have an id.
      $read_file_id = $read_file->dbID;
    }
    
    my $experiment = $current_object->get_Experiment;
    my $experiment_id = undef;
    
    if (! defined $experiment) {
      throw("Experiment must be set!")
    }
    
    $experiment_id = $experiment->dbID;
    
    $sth_store_object->bind_param( 1, $read_file_id,                         SQL_INTEGER);
    $sth_store_object->bind_param( 2, $experiment_id,                        SQL_INTEGER);
    $sth_store_object->bind_param( 3, $current_object->biological_replicate, SQL_INTEGER);
    $sth_store_object->bind_param( 4, $current_object->technical_replicate,  SQL_INTEGER);
    
    $sth_store_object->execute;
    $current_object->dbID( $self->last_insert_id );
  }
  return;
}

1;
