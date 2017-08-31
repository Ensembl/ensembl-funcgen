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

package Bio::EnsEMBL::Funcgen::DBSQL::ReadFileAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use DBI qw(:sql_types);

use base 'Bio::EnsEMBL::DBSQL::BaseAdaptor';

sub _tables {
  return (
    ['read_file', 'rf'  ],
  );
}

sub _columns {
  my $self = shift;
  
  return qw(
      rf.read_file_id
      rf.name
      rf.analysis_id
      rf.is_paired_end
      rf.paired_with
      rf.file_size
      rf.read_length
  );
}

sub _default_where_clause {
  return '';
}

sub fetch_by_name {
  my $self = shift;
  my $name = shift;

  my $constraint = "rf.name = ?";
  $self->bind_param_generic_fetch($name, SQL_VARCHAR);
  
  my $object_list = $self->generic_fetch($constraint);
  
  if (!$object_list || @$object_list == 0) {
    return;
  }
  if (@$object_list != 1) {
    throw("Found ". @$object_list ." read files with the same name!");
  }
  return $object_list->[0];
}

sub _objs_from_sth {
  my ($self, $sth) = @_;

  my(
    $sth_fetched_dbID,
    $sth_fetched_name,
    $sth_fetched_is_paired_end,
    $sth_fetched_paired_with,
    $sth_fetched_file_size,
    $sth_fetched_read_length,
    $sth_fetched_analysis_id,
  );
  
  $sth->bind_columns (
    \$sth_fetched_dbID,
    \$sth_fetched_name,
    \$sth_fetched_analysis_id,
    \$sth_fetched_is_paired_end,
    \$sth_fetched_paired_with,
    \$sth_fetched_file_size,
    \$sth_fetched_read_length,
  );
  
  use Bio::EnsEMBL::Funcgen::ReadFile;
  
  my $analysis_adaptor                             = $self->db->get_AnalysisAdaptor;
  my $read_file_experimental_configuration_adaptor = $self->db->get_ReadFileExperimentalConfigurationAdaptor;
  
  my @return_object_list;
  
  ROW: while ( $sth->fetch() ) {
  
    my $analysis = $analysis_adaptor->fetch_by_dbID($sth_fetched_analysis_id);

    my $current_object = Bio::EnsEMBL::Funcgen::ReadFile->new(
      -db                                   => $self->db,
      -dbID                                 => $sth_fetched_dbID,
      -name                                 => $sth_fetched_name,
      -is_paired_end                        => $sth_fetched_is_paired_end,
      -paired_with                          => $sth_fetched_paired_with,
      -file_size                            => $sth_fetched_file_size,
      -read_length                          => $sth_fetched_read_length,
      -analysis                             => $analysis,

    );
    push @return_object_list, $current_object;
  }
  return \@return_object_list;
}

sub store {
  my ($self, @object) = @_;
  
  my $sth_store_object = $self->prepare("
    INSERT INTO read_file (
      name,
      analysis_id,
      is_paired_end,
      paired_with,
      file_size,
      read_length
    ) VALUES (?, ?, ?, ?, ?, ?)"
  );
  
  my $read_file_experimental_configuration_adaptor = $self->db->get_ReadFileExperimentalConfigurationAdaptor;
  
  foreach my $current_object (@object) {
  
    if (! defined $current_object) {
      throw("Got undefined object!");
    }
  
    $sth_store_object->bind_param( 1, $current_object->name,                    SQL_VARCHAR);
    $sth_store_object->bind_param( 2, $current_object->get_Analysis->dbID,      SQL_INTEGER);
    $sth_store_object->bind_param( 3, $current_object->is_paired_end,           SQL_INTEGER);
    $sth_store_object->bind_param( 4, $current_object->paired_with,             SQL_INTEGER);
    $sth_store_object->bind_param( 5, $current_object->file_size,               SQL_INTEGER);
    $sth_store_object->bind_param( 6, $current_object->read_length,             SQL_INTEGER);
    
    $sth_store_object->execute;
    $current_object->dbID( $self->last_insert_id );
  }
  return;
}


1;
