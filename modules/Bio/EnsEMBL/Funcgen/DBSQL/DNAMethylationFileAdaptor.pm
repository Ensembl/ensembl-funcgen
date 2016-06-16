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
    ['analysis_description',  'ad' ],
    ['feature_type',          'ft' ],
    ['dbfile_registry',       'dr' ],
  );
}

sub _columns {
  my $self = shift;
  
  return qw(
    eff.external_feature_file_id
    eff.name
    a.logic_name
    ad.description
    ad.display_label
    ft.so_accession
    ft.so_name
    ft.class
    ft.name
    dr.path
    dr.file_type
  );
}

sub _default_where_clause {
  return 'eff.analysis_id = a.analysis_id'
    . ' and eff.analysis_id = ad.analysis_id'
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
    $sth_fetched_a_logic_name,
    $sth_fetched_ad_description,
    $sth_fetched_ad_display_label,
    $sth_fetched_ft_so_accession,
    $sth_fetched_ft_so_name,
    $sth_fetched_ft_class,
    $sth_fetched_ft_name,
    $sth_fetched_dr_path,
    $sth_fetched_dr_file_type,
  );

  $sth->bind_columns (
    \$sth_fetched_dbID,
    \$sth_fetched_eff_name,
    \$sth_fetched_a_logic_name,
    \$sth_fetched_ad_description,
    \$sth_fetched_ad_display_label,
    \$sth_fetched_ft_so_accession,
    \$sth_fetched_ft_so_name,
    \$sth_fetched_ft_class,
    \$sth_fetched_ft_name,
    \$sth_fetched_dr_path,
    \$sth_fetched_dr_file_type,
  );
  
  use Bio::EnsEMBL::Funcgen::DNAMethylationFile;
  
  my @return_objects;
  ROW: while ( $sth->fetch() ) {
    my $crispr_sites_file = Bio::EnsEMBL::Funcgen::DNAMethylationFile->new(
      -dbID          => $sth_fetched_dbID,
      -name          => $sth_fetched_eff_name,
      -logic_name    => $sth_fetched_a_logic_name,
      -description   => $sth_fetched_ad_description,
      -display_label => $sth_fetched_ad_display_label,
      -so_accession  => $sth_fetched_ft_so_accession,
      -so_name       => $sth_fetched_ft_so_name,
      -feature_class => $sth_fetched_ft_class,
      -feature_name  => $sth_fetched_ft_name,
      -file          => $sth_fetched_dr_path,
      -file_type     => $sth_fetched_dr_file_type,
    );
    push @return_objects, $crispr_sites_file
  }
  return \@return_objects;
}

1;
