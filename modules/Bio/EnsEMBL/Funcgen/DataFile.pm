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

=head1 SYNOPSIS

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::DataFile;

use strict;

use base 'Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality';

sub _constructor_parameters {
  return {
    data_file_id => 'data_file_id',
    table_id     => 'table_id',
    table_name   => 'table_name',
    path         => 'path',
    file_type    => 'file_type',
    md5sum       => 'md5sum',
  };
}

sub _simple_accessors {
  return [
    { method_name => 'data_file_id', hash_key => 'data_file_id' },
    { method_name => 'table_id',     hash_key => 'table_id'     },
    { method_name => 'table_name',   hash_key => 'table_name'   },
    { method_name => 'path',         hash_key => 'path'         },
    { method_name => 'file_type',    hash_key => 'file_type'    },
    { method_name => 'md5sum',       hash_key => 'md5sum'       },
  ]
}

1;
