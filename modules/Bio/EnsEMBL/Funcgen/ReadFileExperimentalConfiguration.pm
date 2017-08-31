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

package Bio::EnsEMBL::Funcgen::ReadFileExperimentalConfiguration;

use strict;
use warnings;

use base 'Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality';

sub _constructor_parameters {
  return {
    technical_replicate  => 'technical_replicate',
    biological_replicate => 'biological_replicate',

    read_file  => 'set_ReadFile',
    experiment => 'set_Experiment',
  }
}

sub _simple_accessors {
  return [
    { method_name => 'biological_replicate', hash_key => '_biological_replicate', },
    { method_name => 'technical_replicate',  hash_key => '_technical_replicate',  },
  ]
}

sub _get_methods {
  return [
    { method_name => 'get_Experiment',      hash_key    => 'experiment',    },
    { method_name => 'get_ReadFile',        hash_key    => 'read_file',     },
  ]
}

sub _set_methods {
  return [
    {
      method_name   => 'set_Experiment',
      expected_type => 'Bio::EnsEMBL::Funcgen::Experiment',
      hash_key      => 'experiment',
    },
    {
      method_name   => 'set_ReadFile',
      expected_type => 'Bio::EnsEMBL::Funcgen::ReadFile',
      hash_key      => 'read_file'
    },
  ]
}

1;

