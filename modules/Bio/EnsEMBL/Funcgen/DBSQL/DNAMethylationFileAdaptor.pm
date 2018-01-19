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

Bio::EnsEMBL::DBSQL::Funcgen::DNAMethylationFileAdaptor

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::DNAMethylationFileAdaptor;

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
    ['feature_type',          'ft' ],
    ['data_file',             'dr' ],
  );
}

sub _columns {
  my $self = shift;
  
  return qw(
    eff.external_feature_file_id
    eff.name
    ft.name
    dr.path
    dr.file_type
    a.analysis_id
    ft.feature_type_id
  );
}

sub _default_where_clause {
  return 'eff.analysis_id = a.analysis_id'
    . ' and eff.feature_type_id = ft.feature_type_id'
    . ' and dr.table_name="external_feature_file" and dr.table_id=external_feature_file_id'
    . ' and ft.name = "5mC"'
    ;
}

sub fetch_by_name {
  my $self  = shift;
  my $name  = shift;

  my $constraint = "eff.name = ?";
  
  $self->bind_param_generic_fetch($name,  SQL_VARCHAR);
  my $dna_methylation_file = $self->generic_fetch($constraint);
  
  if (!$dna_methylation_file || @$dna_methylation_file==0) {
    return;
  }
  if (@$dna_methylation_file!=1) {
    throw("Found ". @$dna_methylation_file ." dna methylation files with the same name!");
  }
  return $dna_methylation_file->[0];
}

sub _objs_from_sth {
  my ($self, $sth) = @_;

  my (
    $sth_fetched_dbID,
    $sth_fetched_eff_name,
    $sth_fetched_ft_name,
    $sth_fetched_dr_path,
    $sth_fetched_dr_file_type,
    $sth_fetched_a_analysis_id,
    $sth_fetched_ft_feature_type_id
  );

  $sth->bind_columns (
    \$sth_fetched_dbID,
    \$sth_fetched_eff_name,
    \$sth_fetched_ft_name,
    \$sth_fetched_dr_path,
    \$sth_fetched_dr_file_type,
    \$sth_fetched_a_analysis_id,
    \$sth_fetched_ft_feature_type_id
  );
  
  use Bio::EnsEMBL::Funcgen::DNAMethylationFile;
  
  my $analysis_adaptor = $self->db->get_AnalysisAdaptor();
  my $feature_type_adaptor = $self->db->get_FeatureTypeAdaptor();
  
  my @return_objects;
  ROW: while ( $sth->fetch() ) {
    my $dna_methylation_file = Bio::EnsEMBL::Funcgen::DNAMethylationFile->new(
      -dbID      => $sth_fetched_dbID,
      -name      => $sth_fetched_eff_name,
      -file      => $sth_fetched_dr_path,
      -file_type => $sth_fetched_dr_file_type,
    );
    
    my $analysis = $analysis_adaptor->fetch_by_dbID($sth_fetched_a_analysis_id);
    $dna_methylation_file->_analysis($analysis);

    my $feature_type = $feature_type_adaptor->fetch_by_dbID($sth_fetched_ft_feature_type_id);
    $dna_methylation_file->_feature_type($feature_type);

    push @return_objects, $dna_methylation_file
  }
  return \@return_objects;
}

1;
