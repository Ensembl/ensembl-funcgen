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

package Bio::EnsEMBL::Funcgen::DBSQL::ProbeMappingStatisticAdaptor;

use strict;
use base 'Bio::EnsEMBL::Funcgen::DBSQL::GenericAdaptor';

sub object_class {
    return 'Bio::EnsEMBL::Funcgen::ProbeMappingStatistic';
}

sub _tables {
  return ['probe_mapping_statistic', 'pms']
}

sub fetch_by_statistic {
  my $self      = shift;
  my $statistic = shift;
  
  my $probe_mapping_adaptor = $self->db->get_RegulatoryBuildAdaptor;
  my $probe_mapping = $probe_mapping_adaptor->fetch_current_probe_mapping;
  
  return $self->fetch_single_object(
    'probe_mapping_id = ? and statistic = ?', 
    [ 
      $probe_mapping->dbID,
      $statistic,
    ]
  );
}

sub fetch_average_length_promoter_flanking_region {
  my $self = shift;
  return $self->fetch_by_statistic('average_length_promoter_flanking_region');
}

1;
