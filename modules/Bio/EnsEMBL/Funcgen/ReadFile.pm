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

package Bio::EnsEMBL::Funcgen::ReadFile;

use strict;
use warnings;

use base 'Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality';

sub _constructor_parameters {
  return {
    name           => 'name',
    is_paired_end  => 'is_paired_end',
    paired_with    => 'paired_with',
    file_size      => 'file_size',
    read_length    => 'read_length',
    md5sum         => 'md5sum',
    file           => 'file',
    notes          => 'notes',
    analysis                             => 'set_Analysis',
    read_file_experimental_configuration => 'set_ReadFileExperimentalConfiguration',
  };
}

sub _simple_accessors {
  return [
    { method_name => 'name',           hash_key => '_name',           },
    { method_name => 'is_paired_end',  hash_key => '_is_paired_end',  },
    { method_name => 'paired_with',    hash_key => '_paired_with',    },
    { method_name => 'file_size',      hash_key => '_file_size',      },
    { method_name => 'read_length',    hash_key => '_read_length',    },
    { method_name => 'md5sum',         hash_key => '_md5sum',         },
    { method_name => 'file',           hash_key => '_file',           },
    { method_name => 'notes',          hash_key => '_notes',          },
  ]
}

sub _get_methods {
  return [
    {
      method_name => 'get_Analysis',
      hash_key    => 'analysis',
    },
    {
      method_name => 'get_ReadFileExperimentalConfiguration',
      hash_key    => 'read_file_experimental_configuration',
    },
  ]
}

sub _set_methods {
  return [
    {
      method_name   => 'set_Analysis',
      expected_type => 'Bio::EnsEMBL::Analysis',
      hash_key      => 'analysis',
    },
    {
      method_name   => 'set_ReadFileExperimentalConfiguration',
      expected_type => 'Bio::EnsEMBL::Funcgen::ReadFileExperimentalConfiguration',
      hash_key      => 'read_file_experimental_configuration',
    },
  ]
}

1;
