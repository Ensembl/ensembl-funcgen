=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

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

=head1 SYNOPSIS

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::RegulatoryBuildStatisticsAdaptor;

use strict;
use base 'Bio::EnsEMBL::Funcgen::DBSQL::GenericAdaptor';

sub object_class {
    return 'Bio::EnsEMBL::Funcgen::RegulatoryBuildStatistic';
}

sub _tables {
  return ['regulatory_build_statistics', 'rbs']
}

sub fetch_by_statistic {
  my $self      = shift;
  my $statistic = shift;
  
  my $regulatory_build_adaptor = $self->db->get_RegulatoryBuildAdaptor;
  my $regulatory_build = $regulatory_build_adaptor->fetch_current_regulatory_build;
  
  return $self->fetch_single_object(
    'regulatory_build_id = ? and statistic = ?', 
    [ 
      $regulatory_build->dbID,
      $statistic,
    ]
  );
}

sub fetch_average_length_promoter_flanking_region {
  my $self = shift;
  return $self->fetch_by_statistic('average_length_promoter_flanking_region');
}



sub fetch_sum_length_promoter {
  my $self = shift;
  return $self->fetch_by_statistic('sum_length_promoter');
}

sub fetch_sum_length_enhancer {
  my $self = shift;
  return $self->fetch_by_statistic('sum_length_enhancer');
}

sub fetch_sum_length_promoter_flanking_region {
  my $self = shift;
  return $self->fetch_by_statistic('sum_length_promoter_flanking_region');
}

sub fetch_sum_length_transcription_factor_binding_site {
  my $self = shift;
  return $self->fetch_by_statistic('sum_length_transcription_factor_binding_site');
}

sub fetch_sum_length_open_chromatin {
  my $self = shift;
  return $self->fetch_by_statistic('sum_length_open_chromatin');
}

sub fetch_sum_length_ctcf_binding_site {
  my $self = shift;
  return $self->fetch_by_statistic('sum_length_ctcf_binding_site');
}



sub fetch_number_promoter {
  my $self = shift;
  return $self->fetch_by_statistic('number_promoter');
}

sub fetch_number_enhancer {
  my $self = shift;
  return $self->fetch_by_statistic('number_enhancer');
}

sub fetch_number_promoter_flanking_region {
  my $self = shift;
  return $self->fetch_by_statistic('number_promoter_flanking_region');
}

sub fetch_number_transcription_factor_binding_site {
  my $self = shift;
  return $self->fetch_by_statistic('number_transcription_factor_binding_site');
}

sub fetch_number_open_chromatin {
  my $self = shift;
  return $self->fetch_by_statistic('number_open_chromatin');
}

sub fetch_number_ctcf_binding_site {
  my $self = shift;
  return $self->fetch_by_statistic('number_ctcf_binding_site');
}




sub fetch_average_length_promoter {
  my $self = shift;
  return $self->fetch_by_statistic('average_length_promoter');
}

sub fetch_average_length_enhancer {
  my $self = shift;
  return $self->fetch_by_statistic('average_length_enhancer');
}

sub fetch_average_length_promoter_flanking_region {
  my $self = shift;
  return $self->fetch_by_statistic('average_length_promoter_flanking_region');
}

sub fetch_average_length_transcription_factor_binding_site {
  my $self = shift;
  return $self->fetch_by_statistic('average_length_transcription_factor_binding_site');
}

sub fetch_average_length_open_chromatin {
  my $self = shift;
  return $self->fetch_by_statistic('average_length_open_chromatin');
}

sub fetch_average_length_ctcf_binding_site {
  my $self = shift;
  return $self->fetch_by_statistic('average_length_ctcf_binding_site');
}


1;
