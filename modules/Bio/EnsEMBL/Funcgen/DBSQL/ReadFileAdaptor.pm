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

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::ReadFileAdaptor;

use strict;
use Bio::EnsEMBL::Utils::Exception qw( throw );
use base 'Bio::EnsEMBL::Funcgen::DBSQL::GenericAdaptor';

sub object_class {
  return 'Bio::EnsEMBL::Funcgen::ReadFile';
}

sub _tables {
  return (
    ['read_file', 'rf'  ],
  );
}

sub fetch_all_by_Experiment {
  my $self = shift;
  my $experiment = shift;
  
  my $read_file_experimental_configuration_adaptor
    = $self->db->get_ReadFileExperimentalConfigurationAdaptor;
  
  my $read_file_experimental_configuration_list
    = $read_file_experimental_configuration_adaptor
      ->fetch_all_by_Experiment($experiment);

  my @all_read_files_for_experiment;
  
  foreach my $read_file_experimental_configuration 
    (@$read_file_experimental_configuration_list) {
      
      my $read_file 
        = $self->fetch_by_ReadFileExperimentalConfiguration
          ($read_file_experimental_configuration);
      
      push @all_read_files_for_experiment, $read_file;
  }
  return \@all_read_files_for_experiment;
}

# Here, the read_file_experimental_configuration is used like a filter in 
# the database.
#
# Find the read file in the database that is described by this.
#
# That is why the get_ReadFile can't be assumed to be populated and is
# not used.
#
sub fetch_by_ReadFileExperimentalConfiguration {
  my $self = shift;
  my $read_file_experimental_configuration = shift;
  
  if (! defined $read_file_experimental_configuration) {
    throw("Read file experimental configuration parameter was undefined!");
  }
  
#   my $read_file = $read_file_experimental_configuration->get_ReadFile;
#   if (defined $read_file) {
#     return $read_file;
#   }
  
  my $read_file_experimental_configuration_adaptor
    = $self->db->get_ReadFileExperimentalConfigurationAdaptor;

  my $read_file_experimental_configuration_from_db 
    = $read_file_experimental_configuration_adaptor
      ->fetch_by_ReadFileExperimentalConfiguration(
        $read_file_experimental_configuration
      );
  
  if (! defined $read_file_experimental_configuration_from_db) {
    use Data::Dumper;
    throw(
      "Couldn't find read file experimental configuration in database:\n\n"
      . Dumper($read_file_experimental_configuration)
    );
  }
  
  return $read_file_experimental_configuration_from_db->get_ReadFile;
}

sub _load_dependencies {
    my $self = shift;
    my $read_file = shift;
    
    my $analysis_id = $read_file->_analysis_id;
    
    my $analysis_adaptor = $self->db->get_AnalysisAdaptor;
    my $analysis = $analysis_adaptor->fetch_by_dbID($analysis_id);
    
    $read_file->set_Analysis($analysis);

    return;
}

sub _store_dependencies {
  my $self = shift;
  my $read_file = shift;
  
  my $analysis = $read_file->get_Analysis;
  my $analysis_id = $analysis->dbID;
  
  if (! defined $analysis_id) {
    throw("Analysis must exist in the database already!");
  }

  my $read_file_id = $read_file->dbID;
  
  $self->sql_helper->execute_update(
    -SQL      => '
      update 
        read_file
      set 
        analysis_id = ? 
      where 
        read_file_id = ?
    ',
    -PARAMS => [ $analysis_id, $read_file_id ],
  );
  return;
}

1;
