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

  Bio::EnsEMBL::Funcgen::DBSQL::ProbeFeatureTranscriptMappingAdaptor

=head1 SYNOPSIS

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::ProbeFeatureTranscriptMappingAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use DBI qw(:sql_types);

use vars '@ISA';
@ISA    = qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

sub _tables {
  return ['probe_feature_transcript', 'pt'];
}

sub _columns {
  my $self = shift;
  
  return qw(
    pt.probe_feature_transcript_id
    pt.probe_feature_id
    pt.stable_id
    pt.description
  );
}

sub fetch_all_by_transcript_stable_id {
  my $self  = shift;
  my $transcript_stable_id = shift;

  my $constraint = "pt.stable_id = ?";
  
  $self->bind_param_generic_fetch($transcript_stable_id, SQL_VARCHAR);
  my $mapping = $self->generic_fetch($constraint);
  
  return $mapping;
}

sub fetch_all_by_probe_feature_id {
  my $self  = shift;
  my $probe_feature_id = shift;

  my $constraint = "pt.probe_feature_id = ?";
  
  $self->bind_param_generic_fetch($probe_feature_id, SQL_VARCHAR);
  my $mapping = $self->generic_fetch($constraint);
  
  return $mapping;
}

sub _objs_from_sth {
  my ($self, $sth) = @_;

  my (
    $sth_fetched_dbID,
    $sth_fetched_probe_feature_id,
    $sth_fetched_stable_id,
    $sth_fetched_description,
  );

  $sth->bind_columns (
    \$sth_fetched_dbID,
    \$sth_fetched_probe_feature_id,
    \$sth_fetched_stable_id,
    \$sth_fetched_description,
  );
  
  use Bio::EnsEMBL::Funcgen::ProbeTranscriptMapping;

  my @return_objects;
  ROW: while ( $sth->fetch() ) {
  
    my $probe_feature_transcript_mapping = Bio::EnsEMBL::Funcgen::ProbeFeatureTranscriptMapping->new(
      -dbID             => $sth_fetched_dbID,
      -probe_feature_id => $sth_fetched_probe_feature_id,
      -stable_id        => $sth_fetched_stable_id,
      -description      => $sth_fetched_description,
      -adaptor          => $self->db,
    );
    push @return_objects, $probe_feature_transcript_mapping
  }
  return \@return_objects;
}

1;
