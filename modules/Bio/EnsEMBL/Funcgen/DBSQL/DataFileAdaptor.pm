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

package Bio::EnsEMBL::Funcgen::DBSQL::DataFileAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw  );
use base 'Bio::EnsEMBL::Funcgen::DBSQL::GenericAdaptor';

sub object_class {
    return 'Bio::EnsEMBL::Funcgen::DataFile';
}

sub _tables {
  return ['data_file', 'df']
}

sub fetch_by_path {
  my $self = shift;
  my $path = shift;
  
  my $data_file_objects = $self->fetch_all('path = ?', [ $path ]);
  
  if (@$data_file_objects == 0) {
    return;
  }
  if (@$data_file_objects > 1) {
    throw("Found more than one data file with path $path!");
  }
  return $data_file_objects->[0];
}

sub _fetch_all_by_table_name {
    my ($self, $table_name) = @_;

    my %valid_table_names = (
        'motif_feature'         => 1,
        'alignment'             => 1,
        'segmentation_file'     => 1,
        'external_feature_file' => 1
    );

    if (!exists $valid_table_names{$table_name}) {
        throw('Please provide one of the following valid table_name parameters: '
            . join(',', keys %valid_table_names));
    }
    my $rows = $self->fetch_all('table_name = ?', [ $table_name ]);

    if (scalar @{$rows} > 0){
        return $rows;
    }
    else {
        return [];
    }
}

1;
