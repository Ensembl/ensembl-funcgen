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


=head1 NAME

Bio::EnsEMBL::Funcgen::DBSQL::ReadFileExperimentalConfigurationAdaptor

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::ReadFileExperimentalConfigurationAdaptor;

use strict;
use base 'Bio::EnsEMBL::Funcgen::DBSQL::GenericAdaptor';

sub object_class {
  return 'Bio::EnsEMBL::Funcgen::ReadFileExperimentalConfiguration';
}

sub _tables {
  return (
    ['read_file_experimental_configuration', 'rfec'  ],
  );
}

sub _load_dependencies {
    my $self = shift;
    my $read_file_experimental_configuration = shift;
    my $read_file_id = $read_file_experimental_configuration->_read_file_id;
    
    my $read_file_adaptor = $self->db->get_ReadFileAdaptor;
    my $read_file = $read_file_adaptor->fetch_by_dbID($read_file_id);
    
    $read_file_experimental_configuration->set_ReadFile($read_file);
    return;
}

sub _store_dependencies {

  my $self = shift;
  my $read_file_experimental_configuration = shift;
  
  my $read_file = $read_file_experimental_configuration->get_ReadFile;
  my $read_file_id = $read_file->dbID;
  my $read_file_experimental_configuration_id = $read_file_experimental_configuration->dbID;
  
  if (! defined $read_file_id) {
    my $read_file_adaptor = $self->db->get_ReadFileAdaptor;
    $read_file_adaptor->store($read_file);
    $read_file_id = $read_file->dbID;
  }

  $self->sql_helper->execute_update(
    -SQL      => '
      update 
        read_file_experimental_configuration 
      set 
        read_file_id = ? 
      where 
          read_file_experimental_configuration_id = ?
    ',
    -PARAMS => [ $read_file_id, $read_file_experimental_configuration_id ],
  );
  return;
}

sub fetch_all_by_read_file_id {

  my $self = shift;
  my $read_file_id = shift;
  
  my $reads = $self->fetch_all(
    "read_file_id = " . $read_file_id
  );
  return $reads;
}

sub fetch_all_by_experiment_id {

  my $self          = shift;
  my $experiment_id = shift;
  
  my $read_file_experimental_configuration_list = $self->fetch_all(
    "experiment_id = " . $experiment_id
  );
  return $read_file_experimental_configuration_list;
}

sub fetch_all_by_Experiment {
  my $self = shift;
  my $experiment = shift;
  my $experiment_id = $experiment->dbID;
  return $self->fetch_all_by_experiment_id($experiment_id);
}

sub fetch_all_technical_replicates_by_Experiment_and_biological_replicate_number {
  my $self = shift;
  my $experiment                  = shift;
  my $biological_replicate_number = shift;
  
  my $experiment_id = $experiment->dbID;
  
  my $read_file_experimental_configuration_list = $self->fetch_all(
    "experiment_id = $experiment_id"
    . " and biological_replicate = $biological_replicate_number"
  );
  
  return $read_file_experimental_configuration_list;
}

=head2 fetch_all_biological_replicate_numbers_from_Experiment

  Description: Convenience method
  Returntype : ArrayRef[Int]
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
sub fetch_all_biological_replicate_numbers_from_Experiment {

  my $self       = shift;
  my $experiment = shift;
  
  my $experiment_id = $experiment->dbID;

  my @biological_replicate_numbers;

  $self->sql_helper->execute_no_return(
    -SQL          => '
      select 
        distinct biological_replicate 
      from 
        read_file_experimental_configuration 
      where 
        experiment_id = ? 
      order by 
        biological_replicate
    ',
    -PARAMS       => [ $experiment_id ],
    -USE_HASHREFS => 1,
    -CALLBACK     => sub {
        my $row = shift;
        my $biological_replicate = $row->{biological_replicate};
        push @biological_replicate_numbers, $biological_replicate;
        return;
      },
  );
  return \@biological_replicate_numbers;
}

sub count_biological_replicates_from_Experiment {

    my $self = shift;
    my $experiment = shift;
    
    my $biological_replicate_to_counts = $self->count_all(
        "experiment_id = " . $experiment->dbID,
        [ 'biological_replicate' ]
    );
    my @biological_replicates = keys %$biological_replicate_to_counts;
    return scalar @biological_replicates;
}

sub count_technical_replicates_from_Experiment {

    my $self = shift;
    my $experiment = shift;
    
    my $technical_replicate_to_counts = $self->count_all(
        "experiment_id = " . $experiment->dbID,
        [ 'technical_replicate' ]
    );
    my @technical_replicates = keys %$technical_replicate_to_counts;
    return scalar @technical_replicates;
}

1;
