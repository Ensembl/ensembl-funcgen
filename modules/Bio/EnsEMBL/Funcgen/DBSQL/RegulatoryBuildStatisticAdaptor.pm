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

=head1 SYNOPSIS

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::RegulatoryBuildStatisticAdaptor;

use strict;
use base 'Bio::EnsEMBL::Funcgen::DBSQL::GenericAdaptor';
use Bio::EnsEMBL::Utils::Exception qw( throw warning );

sub object_class {
    return 'Bio::EnsEMBL::Funcgen::RegulatoryBuildStatistic';
}

sub _tables {
  return ['regulatory_build_statistic', 'rbs']
}

sub fetch_by_statistic {
  my $self      = shift;
  my $statistic = shift;

  my $error_result = Bio::EnsEMBL::Funcgen::RegulatoryBuildStatistic->new(
        statistic => $statistic,
        value     => -1,
  );
  my $result = $self->_generic_fetch_by_statistic($statistic);
  
  if (! defined $result->[0]) {
    warning("Statistic $statistic does not exist!");
    $result->[0] = $error_result;
  }
  if (! defined $result->[0]->value) {
    warning("Statistic is not defined!");
    $result->[0] = $error_result;
  }
  if (@$result > 1) {
    warning("Statistic $statistic is not unique!");
    $result->[0] = $error_result;
  }
  return $result->[0];
}

sub _statistic_exists {
  my $self      = shift;
  my $statistic = shift;
  
  my $result = $self->_generic_fetch_by_statistic($statistic);
  my $statistic_exists = @$result != 0;
  return $statistic_exists
}

sub _generic_fetch_by_statistic {
  my $self      = shift;
  my $statistic = shift;
  
  my $regulatory_build_adaptor = $self->db->get_RegulatoryBuildAdaptor;
  my $regulatory_build = $regulatory_build_adaptor->fetch_current_regulatory_build;
  
   return $self->fetch_all(
      'regulatory_build_id = ? and statistic = ?', 
      [ 
        $regulatory_build->dbID,
        $statistic,
      ]
    )
  ;
}

sub fetch_num_epigenomes_in_regulatory_build {

  my $self = shift;

  my $regulatory_build_adaptor = $self->db->get_RegulatoryBuildAdaptor;
  my $regulatory_build = $regulatory_build_adaptor->fetch_current_regulatory_build;
  my $epigenomes_in_regulatory_build = $regulatory_build->get_all_Epigenomes;
  
  my $num_epigenomes_in_regulatory_build = scalar @$epigenomes_in_regulatory_build;
  
  use Bio::EnsEMBL::Funcgen::RegulatoryBuildStatistic;
  return Bio::EnsEMBL::Funcgen::RegulatoryBuildStatistic->new(
    -value => $num_epigenomes_in_regulatory_build,
  );
}

sub fetch_regulatory_build_overlap_percent {
  my $self = shift;
  return $self->fetch_by_statistic('regulatory_build_overlap_percent');
}
sub fetch_ctcf_overlap_percent {
  my $self = shift;
  return $self->fetch_by_statistic('ctcf_overlap_percent');
}
sub fetch_enhancer_overlap_percent {
  my $self = shift;
  return $self->fetch_by_statistic('enhancer_overlap_percent');
}
sub fetch_tf_binding_overlap_percent {
  my $self = shift;
  return $self->fetch_by_statistic('tf_binding_overlap_percent');
}
sub fetch_promoter_overlap_percent {
  my $self = shift;
  return $self->fetch_by_statistic('promoter_overlap_percent');
}
sub fetch_promoter_flanking_overlap_percent {
  my $self = shift;
  return $self->fetch_by_statistic('promoter_flanking_overlap_percent');
}
sub fetch_open_chromatin_overlap_percent {
  my $self = shift;
  return $self->fetch_by_statistic('open_chromatin_overlap_percent');
}

sub fetch_regulatory_build_overlap_bp {
  my $self = shift;
  return $self->fetch_by_statistic('regulatory_build_overlap_bp');
}
sub fetch_ctcf_overlap_bp {
  my $self = shift;
  return $self->fetch_by_statistic('ctcf_overlap_bp');
}
sub fetch_enhancer_overlap_bp {
  my $self = shift;
  return $self->fetch_by_statistic('enhancer_overlap_bp');
}
sub fetch_tf_binding_overlap_bp {
  my $self = shift;
  return $self->fetch_by_statistic('tf_binding_overlap_bp');
}
sub fetch_promoter_overlap_bp {
  my $self = shift;
  return $self->fetch_by_statistic('promoter_overlap_bp');
}
sub fetch_promoter_flanking_overlap_bp {
  my $self = shift;
  return $self->fetch_by_statistic('promoter_flanking_overlap_bp');
}
sub fetch_open_chromatin_overlap_bp {
  my $self = shift;
  return $self->fetch_by_statistic('open_chromatin_overlap_bp');
}

sub fetch_promoter_q0 {
  my $self = shift;
  return $self->fetch_by_statistic('promoter_q0');
}
sub fetch_promoter_q1 {
  my $self = shift;
  return $self->fetch_by_statistic('promoter_q1');
}
sub fetch_promoter_q2 {
  my $self = shift;
  return $self->fetch_by_statistic('promoter_q2');
}
sub fetch_promoter_q3 {
  my $self = shift;
  return $self->fetch_by_statistic('promoter_q3');
}
sub fetch_promoter_q4 {
  my $self = shift;
  return $self->fetch_by_statistic('promoter_q4');
}
sub fetch_promoter_skewness {
  my $self = shift;
  return $self->fetch_by_statistic('promoter_skewness');
}
sub fetch_promoter_kurtosis {
  my $self = shift;
  return $self->fetch_by_statistic('promoter_kurtosis');
}

sub fetch_promoter_flanking_q0 {
  my $self = shift;
  return $self->fetch_by_statistic('promoter_flanking_q0');
}
sub fetch_promoter_flanking_q1 {
  my $self = shift;
  return $self->fetch_by_statistic('promoter_flanking_q1');
}
sub fetch_promoter_flanking_q2 {
  my $self = shift;
  return $self->fetch_by_statistic('promoter_flanking_q2');
}
sub fetch_promoter_flanking_q3 {
  my $self = shift;
  return $self->fetch_by_statistic('promoter_flanking_q3');
}
sub fetch_promoter_flanking_q4 {
  my $self = shift;
  return $self->fetch_by_statistic('promoter_flanking_q4');
}
sub fetch_promoter_flanking_skewness {
  my $self = shift;
  return $self->fetch_by_statistic('promoter_flanking_skewness');
}
sub fetch_promoter_flanking_kurtosis {
  my $self = shift;
  return $self->fetch_by_statistic('promoter_flanking_kurtosis');
}



sub fetch_enhancer_q0 {
  my $self = shift;
  return $self->fetch_by_statistic('enhancer_q0');
}
sub fetch_enhancer_q1 {
  my $self = shift;
  return $self->fetch_by_statistic('enhancer_q1');
}
sub fetch_enhancer_q2 {
  my $self = shift;
  return $self->fetch_by_statistic('enhancer_q2');
}
sub fetch_enhancer_q3 {
  my $self = shift;
  return $self->fetch_by_statistic('enhancer_q3');
}
sub fetch_enhancer_q4 {
  my $self = shift;
  return $self->fetch_by_statistic('enhancer_q4');
}
sub fetch_enhancer_skewness {
  my $self = shift;
  return $self->fetch_by_statistic('enhancer_skewness');
}
sub fetch_enhancer_kurtosis {
  my $self = shift;
  return $self->fetch_by_statistic('enhancer_kurtosis');
}

sub fetch_ctcf_q0 {
  my $self = shift;
  return $self->fetch_by_statistic('ctcf_q0');
}
sub fetch_ctcf_q1 {
  my $self = shift;
  return $self->fetch_by_statistic('ctcf_q1');
}
sub fetch_ctcf_q2 {
  my $self = shift;
  return $self->fetch_by_statistic('ctcf_q2');
}
sub fetch_ctcf_q3 {
  my $self = shift;
  return $self->fetch_by_statistic('ctcf_q3');
}
sub fetch_ctcf_q4 {
  my $self = shift;
  return $self->fetch_by_statistic('ctcf_q4');
}
sub fetch_ctcf_skewness {
  my $self = shift;
  return $self->fetch_by_statistic('ctcf_skewness');
}
sub fetch_ctcf_kurtosis {
  my $self = shift;
  return $self->fetch_by_statistic('ctcf_kurtosis');
}

sub fetch_tf_q0 {
  my $self = shift;
  return $self->fetch_by_statistic('tf_q0');
}
sub fetch_tf_q1 {
  my $self = shift;
  return $self->fetch_by_statistic('tf_q1');
}
sub fetch_tf_q2 {
  my $self = shift;
  return $self->fetch_by_statistic('tf_q2');
}
sub fetch_tf_q3 {
  my $self = shift;
  return $self->fetch_by_statistic('tf_q3');
}
sub fetch_tf_q4 {
  my $self = shift;
  return $self->fetch_by_statistic('tf_q4');
}
sub fetch_tf_skewness {
  my $self = shift;
  return $self->fetch_by_statistic('tf_skewness');
}
sub fetch_tf_kurtosis {
  my $self = shift;
  return $self->fetch_by_statistic('tf_kurtosis');
}

sub fetch_open_chromatin_q0 {
  my $self = shift;
  return $self->fetch_by_statistic('open_chromatin_q0');
}
sub fetch_open_chromatin_q1 {
  my $self = shift;
  return $self->fetch_by_statistic('open_chromatin_q1');
}
sub fetch_open_chromatin_q2 {
  my $self = shift;
  return $self->fetch_by_statistic('open_chromatin_q2');
}
sub fetch_open_chromatin_q3 {
  my $self = shift;
  return $self->fetch_by_statistic('open_chromatin_q3');
}
sub fetch_open_chromatin_q4 {
  my $self = shift;
  return $self->fetch_by_statistic('open_chromatin_q4');
}
sub fetch_open_chromatin_skewness {
  my $self = shift;
  return $self->fetch_by_statistic('open_chromatin_skewness');
}
sub fetch_open_chromatin_kurtosis {
  my $self = shift;
  return $self->fetch_by_statistic('open_chromatin_kurtosis');
}

sub fetch_num_enhancers_overlapping_vista {
  my $self = shift;
  return $self->fetch_by_statistic('num_enhancers_overlapping_vista');
}
sub fetch_total_enhancers_checked_vista {
  my $self = shift;
  return $self->fetch_by_statistic('total_enhancers_checked_vista');
}

sub fetch_num_enhancers_overlapping_fantom {
  my $self = shift;
  return $self->fetch_by_statistic('num_enhancers_overlapping_fantom');
}
sub fetch_total_enhancers_checked_fantom {
  my $self = shift;
  return $self->fetch_by_statistic('total_enhancers_checked_fantom');
}

sub has_stable_id_mapping_statistics {
  my $self = shift;
  
  my $stable_id_mapping_statistics_available
    = $self->_statistic_exists('stable_id_mapping_number_regulatory_features');
  
  return $stable_id_mapping_statistics_available;
}
sub fetch_stable_id_mapping_number_regulatory_features {
  my $self = shift;
  return $self->fetch_by_statistic('stable_id_mapping_number_regulatory_features');
}
sub fetch_stable_id_mapping_new_stable_ids {
  my $self = shift;
  return $self->fetch_by_statistic('stable_id_mapping_new_stable_ids');
}
sub fetch_stable_id_mapping_mapped_stable_ids {
  my $self = shift;
  return $self->fetch_by_statistic('stable_id_mapping_mapped_stable_ids');
}

sub fetch_number_regulatory_features {
  my $self = shift;
  return $self->fetch_by_statistic('number_regulatory_features');
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
