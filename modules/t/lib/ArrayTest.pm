=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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

=cut

package ArrayTest;

use strict;
use warnings;

use Test::More;
use Test::Exception;

use parent qw(Bio::EnsEMBL::Funcgen::Test);

sub parameters :Test(setup) {
    my $self = shift;

    my %mandatory_constructor_parameters = (
        -name   => 'dummy_array_name',
        -vendor => 'dummy_vendor_name'
    );

    $self->{mandatory_constructor_parameters} =
        \%mandatory_constructor_parameters;

    my %optional_constructor_parameters = (
        -format                  => 'dummy_format',
        -type                    => 'dummy_type',
        -description             => 'dummy_description',
        -class                   => 'dummy_class',
        -is_probeset_array       => 1,
        -is_linked_array         => 1,
        -has_sense_interrogation => 1
    );

    my %constructor_parameters = (%mandatory_constructor_parameters,
                                  %optional_constructor_parameters);

    $self->{constructor_parameters} = \%constructor_parameters;

    my %parameters = (
    );

    $self->{parameters} = \%parameters;
}

sub define_expected :Test(setup) {
    my $self = shift;

    $self->{expected} = {
        'name'                    => 'HumanWG_6_V2',
        'format'                  => 'EXPRESSION',
        'vendor'                  => 'ILLUMINA',
        'description'             => '',
        'type'                    => 'OLIGO',
        'class'                   => 'ILLUMINA_WG',
        'is_probeset_array'       => 0,
        'is_linked_array'         => 0,
        'has_sense_interrogation' => 1
    };
}

sub dbIDs_to_fetch {return [ 28 ];}

sub getters_setters {
    return ['is_probeset_array', 'is_linked_array', 'has_sense_interrogation',
    'name', 'type', 'format', 'class', 'vendor', 'description'];
}

1;