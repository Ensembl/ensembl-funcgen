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

Bio::EnsEMBL::DBSQL::Funcgen::RegulatoryBuildAdaptor

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::RegulatoryBuildAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use DBI qw(:sql_types);

use vars '@ISA';
@ISA    = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

sub _tables {
  return (
    ['regulatory_build', 'rb'  ],
  );
}

sub _columns {
  my $self = shift;
  
  return qw(
    rb.regulatory_build_id
    rb.name
    rb.version
    rb.initial_release_date
    rb.last_annotation_update
    rb.feature_type_id
    rb.analysis_id
    rb.is_current
    rb.sample_regulatory_feature_id
  );
}

sub _default_where_clause {
  return '';
}

sub fetch_by_name {
  my $self = shift;
  my $name = shift;

  my $constraint = "rb.name = ?";
  $self->bind_param_generic_fetch($name, SQL_VARCHAR);
  
  my $regulatory_build = $self->generic_fetch($constraint);
  
  if (!$regulatory_build || @$regulatory_build==0) {
    return;
  }
  if (@$regulatory_build!=1) {
    throw("Found ". @$regulatory_build ." current regulatory builds!");
  }
  return $regulatory_build->[0];
}

sub fetch_current_regulatory_build {
  my $self  = shift;

  my $constraint = "rb.is_current = 1";
  my $regulatory_build = $self->generic_fetch($constraint);
  
  if (!$regulatory_build || @$regulatory_build==0) {
    return;
  }
  if (@$regulatory_build!=1) {
    throw("Found ". @$regulatory_build ." current regulatory builds!");
  }
  return $regulatory_build->[0];
}

sub _objs_from_sth {
  my ($self, $sth) = @_;

  my (
    $sth_fetched_dbID,
    $sth_fetched_name,
    $sth_fetched_version,
    $sth_fetched_initial_release_date,
    $sth_fetched_last_annotation_update,
    $sth_fetched_feature_type_id,
    $sth_fetched_analysis_id,
    $sth_fetched_is_current,
    $sth_fetched_sample_regulatory_feature_id,
  );

  $sth->bind_columns (
    \$sth_fetched_dbID,
    \$sth_fetched_name,
    \$sth_fetched_version,
    \$sth_fetched_initial_release_date,
    \$sth_fetched_last_annotation_update,
    \$sth_fetched_feature_type_id,
    \$sth_fetched_analysis_id,
    \$sth_fetched_is_current,
    \$sth_fetched_sample_regulatory_feature_id,
  );
  
  use Bio::EnsEMBL::Funcgen::RegulatoryBuild;
  
  my @return_regulatory_build;
  
  ROW: while ( $sth->fetch() ) {

    my $current_regulatory_build = Bio::EnsEMBL::Funcgen::RegulatoryBuild->new(
      -db                           => $self->db,
      -dbID                         => $sth_fetched_dbID,
      -name                         => $sth_fetched_name,
      -version                      => $sth_fetched_version,
      -initial_release_date         => $sth_fetched_initial_release_date,
      -last_annotation_update       => $sth_fetched_last_annotation_update,
      -feature_type_id              => $sth_fetched_feature_type_id,
      -analysis_id                  => $sth_fetched_analysis_id,
      -is_current                   => $sth_fetched_is_current,
      -sample_regulatory_feature_id => $sth_fetched_sample_regulatory_feature_id,
      
    );
    push @return_regulatory_build, $current_regulatory_build;
  }
  return \@return_regulatory_build;
}

sub _unset_all_current {
  my $self = shift;
  my $sth_unset_current = $self->prepare('update regulatory_build set is_current=0');
  $sth_unset_current->execute;
}

sub store {
  my ($self, @regulatory_build) = @_;
  
  my $sth_store_regulatory_build = $self->prepare("
    INSERT INTO regulatory_build (
      name,
      version,
      initial_release_date,
      last_annotation_update,
      feature_type_id,
      analysis_id,
      is_current,
      sample_regulatory_feature_id
    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)"
  );
  
  foreach my $current_regulatory_build (@regulatory_build) {
  
    if ($current_regulatory_build->is_current) {
      $self->_unset_all_current;
    }
    $sth_store_regulatory_build->bind_param( 1, $current_regulatory_build->name,                         SQL_VARCHAR);
    $sth_store_regulatory_build->bind_param( 2, $current_regulatory_build->version,                      SQL_VARCHAR);
    $sth_store_regulatory_build->bind_param( 3, $current_regulatory_build->initial_release_date,         SQL_VARCHAR);
    $sth_store_regulatory_build->bind_param( 4, $current_regulatory_build->last_annotation_update,       SQL_VARCHAR);
    $sth_store_regulatory_build->bind_param( 5, $current_regulatory_build->feature_type_id,              SQL_INTEGER);
    $sth_store_regulatory_build->bind_param( 6, $current_regulatory_build->analysis_id,                  SQL_INTEGER);
    $sth_store_regulatory_build->bind_param( 7, $current_regulatory_build->is_current,                   SQL_TINYINT);
    $sth_store_regulatory_build->bind_param( 8, $current_regulatory_build->sample_regulatory_feature_id, SQL_TINYINT);
    
    $sth_store_regulatory_build->execute;
    $current_regulatory_build->dbID( $self->last_insert_id );
  }
  return;
}

sub update {
  my ($self, $regulatory_build) = @_;
  
  if (! defined $regulatory_build->dbID) {
    throw("The database id must be set!");
  }
  
  my $sth = $self->prepare("
    update regulatory_build set 
      name = ?,
      version = ?,
      initial_release_date = ?,
      last_annotation_update = ?,
      feature_type_id = ?,
      analysis_id = ?,
      is_current = ?
      sample_regulatory_feature_id = ?
    where regulatory_build_id = ?
  ");
  
  if ($regulatory_build->is_current) {
    $self->_unset_all_current;
  }
  
  $sth->bind_param(1, $regulatory_build->name);
  $sth->bind_param(2, $regulatory_build->version);
  $sth->bind_param(3, $regulatory_build->initial_release_date);
  $sth->bind_param(4, $regulatory_build->last_annotation_update);
  $sth->bind_param(5, $regulatory_build->feature_type_id);
  $sth->bind_param(6, $regulatory_build->analysis_id);
  $sth->bind_param(7, $regulatory_build->is_current);
  $sth->bind_param(8, $regulatory_build->sample_regulatory_feature_id);
  $sth->bind_param(9, $regulatory_build->dbID);
  
  $sth->execute;
  
  return;
}

1;
