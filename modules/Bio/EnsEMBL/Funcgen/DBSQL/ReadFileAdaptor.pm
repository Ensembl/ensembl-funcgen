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
