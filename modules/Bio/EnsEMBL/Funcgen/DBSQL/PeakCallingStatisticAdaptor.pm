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

package Bio::EnsEMBL::Funcgen::DBSQL::PeakCallingStatisticAdaptor;

use strict;
use Carp;
use base 'Bio::EnsEMBL::Funcgen::DBSQL::GenericAdaptor';

sub object_class {
  return 'Bio::EnsEMBL::Funcgen::PeakCallingStatistic';
}

sub _tables {
  return ['peak_calling_statistic', 'pcs']
}

sub fetch_by_PeakCalling_total_length {
  my $self         = shift;
  my $peak_calling = shift;
  return $self->_fetch_by_PeakCalling_statistic($peak_calling, 'total_length');
}

sub fetch_by_PeakCalling_num_peaks {
  my $self         = shift;
  my $peak_calling = shift;
  return $self->_fetch_by_PeakCalling_statistic($peak_calling, 'num_peaks');
}

sub fetch_by_PeakCalling_average_length {
  my $self         = shift;
  my $peak_calling = shift;
  return $self->_fetch_by_PeakCalling_statistic($peak_calling, 'average_length');
}

sub fetch_coverage_percent_by_FeatureType {
  my $self         = shift;
  my $feature_type = shift;
  
  return $self->fetch_single_object(
    'peak_calling_id is null and statistic = "coverage_percent" and feature_type_id = ? and epigenome_id is null', 
    [ 
      $feature_type->dbID,
    ]
  );

}

sub _fetch_by_PeakCalling_statistic {

  my $self = shift;
  
  my $peak_calling = shift;
  my $statistic    = shift;
  
  return $self->fetch_single_object(
    'peak_calling_id = ? and statistic = ?', 
    [ 
      $peak_calling->dbID,
      $statistic,
    ]
  );
}

1;
