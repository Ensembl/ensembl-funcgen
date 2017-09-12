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

Bio::EnsEMBL::DBSQL::Funcgen::SegmentationFileAdaptor

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::SegmentationFileAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use DBI qw(:sql_types);

use vars '@ISA';
@ISA    = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

sub _tables {
  return (
    ['segmentation_file',     'sf' ],
    ['analysis',              'a'  ],
    ['data_file',             'da' ],
  );
}

sub _columns {
  my $self = shift;
  
  return qw(
    sf.segmentation_file_id
    sf.name
    sf.epigenome_id
    sf.regulatory_build_id
    a.analysis_id
    da.path
    da.file_type
  );
}

sub _default_where_clause {
  return 'sf.analysis_id = a.analysis_id'
    . ' and da.table_name="segmentation_file" and da.table_id=segmentation_file_id'
    ;
}

sub fetch_by_name {
  my $self  = shift;
  my $name  = shift;

  my $constraint = "sf.name = ?";
  
  $self->bind_param_generic_fetch($name,  SQL_VARCHAR);
  my $segmentation_file = $self->generic_fetch($constraint);
  
  if (!$segmentation_file || @$segmentation_file==0) {
    return;
  }
  if (@$segmentation_file!=1) {
    throw("Found ". @$segmentation_file ." segmentation files with the same name!");
  }
  return $segmentation_file->[0];
}

sub _objs_from_sth {
  my ($self, $sth) = @_;

  my (
    $sth_fetched_dbID,
    $sth_fetched_sf_name,
    $sth_fetched_sf_epigenome_id,
    $sth_fetched_sf_regulatory_build_id,
    $sth_fetched_a_analysis_id,
    $sth_fetched_dr_path,
    $sth_fetched_dr_file_type,
  );

  $sth->bind_columns (
    \$sth_fetched_dbID,
    \$sth_fetched_sf_name,
    \$sth_fetched_sf_epigenome_id,
    \$sth_fetched_sf_regulatory_build_id,
    \$sth_fetched_a_analysis_id,
    \$sth_fetched_dr_path,
    \$sth_fetched_dr_file_type,
  );
  
  use Bio::EnsEMBL::Funcgen::SegmentationFile;
  
  my $analysis_adaptor         = $self->db->get_AnalysisAdaptor();
  my $epigenome_adaptor        = $self->db->get_EpigenomeAdaptor();
  my $regulatory_build_adaptor = $self->db->get_RegulatoryBuildAdaptor();
  
  my @return_objects;
  ROW: while ( $sth->fetch() ) {
    my $segmentation_file = Bio::EnsEMBL::Funcgen::SegmentationFile->new(
      -dbID          => $sth_fetched_dbID,
      -name          => $sth_fetched_sf_name,
      -file          => $sth_fetched_dr_path,
      -file_type     => $sth_fetched_dr_file_type,
    );
    
    my $analysis = $analysis_adaptor->fetch_by_dbID($sth_fetched_a_analysis_id);
    $segmentation_file->_analysis($analysis);

    my $epigenome = $epigenome_adaptor->fetch_by_dbID($sth_fetched_sf_epigenome_id);
    $segmentation_file->_epigenome($epigenome);

    my $regulatory_build = $regulatory_build_adaptor->fetch_by_dbID($sth_fetched_sf_regulatory_build_id);
    $segmentation_file->_regulatory_build($regulatory_build);

    push @return_objects, $segmentation_file
  }
  return \@return_objects;
}

sub store {
  my ($self, @segmentation_file) = @_;

  if (scalar(@segmentation_file) == 0) {
    throw('Must call store with a list of crispr sites file objects');
  }
  foreach my $current_segmentation_file (@segmentation_file) {
    if( ! ref $current_segmentation_file || ! $current_segmentation_file->isa('Bio::EnsEMBL::Funcgen::SegmentationFile') ) {
      throw('Type error');
    }
  }

  my $sth_store_segmentation_file = $self->prepare("
    INSERT ignore INTO segmentation_file (
      name,
      analysis_id,
      epigenome_id,
      regulatory_build_id
    )
    VALUES (?, ?, ?, ?)"
  );
  $sth_store_segmentation_file->{PrintError} = 0;

  my $sth_store_data_file = $self->prepare("
    INSERT ignore INTO data_file (
      table_id,
      table_name,
      path,
      file_type
    )
    VALUES (?, 'segmentation_file', ?, 'BIGBED')"
  );
  
  my $db = $self->db();

  foreach my $current_segmentation_file (@segmentation_file) {

    $current_segmentation_file->adaptor($self);

    $sth_store_segmentation_file->bind_param( 1, $current_segmentation_file->name,                      SQL_VARCHAR);
    $sth_store_segmentation_file->bind_param( 2, $current_segmentation_file->get_Analysis->dbID,        SQL_INTEGER);
    $sth_store_segmentation_file->bind_param( 3, $current_segmentation_file->get_Epigenome->dbID,       SQL_INTEGER);
    $sth_store_segmentation_file->bind_param( 4, $current_segmentation_file->get_RegulatoryBuild->dbID, SQL_INTEGER);
    
    # Store and set dbID
    $sth_store_segmentation_file->execute;
    $current_segmentation_file->dbID( $self->last_insert_id );

    $sth_store_data_file->bind_param( 1, $current_segmentation_file->dbID, SQL_INTEGER);
    $sth_store_data_file->bind_param( 2, $current_segmentation_file->file, SQL_VARCHAR);
    
    $sth_store_data_file->execute;
  }
  return @segmentation_file;
}

1;
