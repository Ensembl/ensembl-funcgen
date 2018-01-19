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

  Bio::EnsEMBL::DBSQL::Funcgen::CrisprSitesFileAdaptor

=head1 SYNOPSIS

  my $crispr_adaptor = Bio::EnsEMBL::Registry->get_adaptor('homo_sapiens', 'funcgen', 'CrisprSitesFile');
  my $crispr_file = $crispr_adaptor->fetch_file;

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::CrisprSitesFileAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use DBI qw(:sql_types);

use vars '@ISA';
@ISA    = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

sub _tables {
  return (
    ['external_feature_file', 'eff'],
    ['analysis',              'a'  ],
    ['data_file',             'da' ],
  );
}

sub _columns {
  my $self = shift;
  
  return qw(
    eff.external_feature_file_id
    eff.name
    da.path
    da.file_type
    a.analysis_id
  );
}

sub _default_where_clause {
  return 'eff.analysis_id = a.analysis_id'
    . ' and da.table_name="external_feature_file" and da.table_id=external_feature_file_id'
    . ' and a.logic_name = "Crispr"'
    ;
}

=head2 fetch_file

  Example    :   my $crispr_adaptor = Bio::EnsEMBL::Registry->get_adaptor('homo_sapiens', 'funcgen', 'CrisprSitesFile');
                 my $crispr_file = $crispr_adaptor->fetch_file;

  Description: Returns the location of the file on the file system. This is a 
               relative path that has to be prefixed with your document root 
               for the species and assembly. E.g.: For human is would be here:
               ftp://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/

  Returntype : String
  Exceptions : none
  Status     : At Risk

=cut

sub fetch_file {
  my $self  = shift;
  my $crispr_sites_file = $self->generic_fetch;
  return $crispr_sites_file->[0];
}

sub _objs_from_sth {
  my ($self, $sth) = @_;

  my (
    $sth_fetched_dbID,
    $sth_fetched_eff_name,
    $sth_fetched_dr_path,
    $sth_fetched_dr_file_type,
    $sth_fetched_a_analysis_id
  );

  $sth->bind_columns (
    \$sth_fetched_dbID,
    \$sth_fetched_eff_name,
    \$sth_fetched_dr_path,
    \$sth_fetched_dr_file_type,
    \$sth_fetched_a_analysis_id
  );
  
  use Bio::EnsEMBL::Funcgen::CrisprSitesFile;
  
  my $analysis_adaptor = $self->db->get_AnalysisAdaptor();
  
  my @return_crispr_file_objects;
  ROW: while ( $sth->fetch() ) {
    my $crispr_sites_file = Bio::EnsEMBL::Funcgen::CrisprSitesFile->new(
      -dbID          => $sth_fetched_dbID,
      -name          => $sth_fetched_eff_name,
      -file          => $sth_fetched_dr_path,
      -file_type     => $sth_fetched_dr_file_type,
    );
    
    my $analysis = $analysis_adaptor->fetch_by_dbID($sth_fetched_a_analysis_id);
    $crispr_sites_file->_analysis($analysis);
    
    push @return_crispr_file_objects, $crispr_sites_file
  }
  if (@return_crispr_file_objects == 0) {
    die("Can't find crispr file in the database!");
  }
  if (@return_crispr_file_objects != 1) {
    die("There should only be one crispr file in the database! Got " . scalar @return_crispr_file_objects . " entries.");
  }
  return \@return_crispr_file_objects;
}

sub store {
  my ($self, @crispr_sites_file) = @_;

  if (scalar(@crispr_sites_file) == 0) {
    throw('Must call store with a list of crispr sites file objects');
  }
  foreach my $current_crispr_sites_file (@crispr_sites_file) {
    if( ! ref $current_crispr_sites_file || ! $current_crispr_sites_file->isa('Bio::EnsEMBL::Funcgen::CrisprSitesFile') ) {
      throw('Feature must be an RegulatoryFeature object');
    }
  }

  my $sth_store_crispr_sites_file = $self->prepare("
    INSERT ignore INTO external_feature_file (
      name,
      analysis_id
    )
    VALUES (?, ?)"
  );

  $sth_store_crispr_sites_file->{PrintError} = 0;

  my $sth_store_data_file = $self->prepare("
    INSERT ignore INTO data_file (
      table_id,
      table_name,
      path,
      file_type
    )
    VALUES (?, 'external_feature_file', ?, 'BIGBED')"
  );
  
  my $db = $self->db();

  foreach my $current_crispr_sites_file (@crispr_sites_file) {

    $current_crispr_sites_file->adaptor($self);

    $sth_store_crispr_sites_file->bind_param( 1, $current_crispr_sites_file->name,            SQL_VARCHAR);
    $sth_store_crispr_sites_file->bind_param( 2, $current_crispr_sites_file->_analysis->dbID, SQL_INTEGER);
    
    # Store and set dbID
    $sth_store_crispr_sites_file->execute;
    $current_crispr_sites_file->dbID( $self->last_insert_id );

    $sth_store_data_file->bind_param( 1, $current_crispr_sites_file->dbID, SQL_INTEGER);
    $sth_store_data_file->bind_param( 2, $current_crispr_sites_file->file, SQL_VARCHAR);
    
    $sth_store_data_file->execute;
  }
  return @crispr_sites_file;
}


1;
