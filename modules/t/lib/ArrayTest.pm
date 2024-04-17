=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2024] EMBL-European Bioinformatics Institute

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

    my $probes = $self->_quick_fetch_all('Probe', [1,2,3]);
    my $probe_sets = $self->_quick_fetch_all('ProbeSet', [1]);
    my $array_chip_ids = [68];
    my $array_chips = $self->_quick_fetch_all('ArrayChip', $array_chip_ids);
    my $design_ids = ['HumanWG_6_V2'];

    $self->{expected} = {
        'name'                    => 'HumanWG_6_V2',
        'format'                  => 'EXPRESSION',
        'vendor'                  => 'ILLUMINA',
        'description'             => '',
        'type'                    => 'OLIGO',
        'class'                   => 'ILLUMINA_WG',
        'is_probeset_array'       => 0,
        'is_linked_array'         => 0,
        'has_sense_interrogation' => 1,
        'probes'                  => $probes,
        'probe_sets'              => $probe_sets,
        'array_chips'             => $array_chips,
        'array_chip_ids'          => $array_chip_ids,
        'design_ids'              => $design_ids,
    };
}

sub dbIDs_to_fetch {return [ 28 ];}

sub getters_setters {
    return ['is_probeset_array', 'is_linked_array', 'has_sense_interrogation',
    'name', 'type', 'format', 'class', 'vendor', 'description'];
}

sub get_all_Probes :Test(1) {
    my $self = shift;

    is_deeply($self->{fetched}->[0]->get_all_Probes,
    $self->{expected}->{probes},
    'get_all_Probes() works');
}

sub get_all_ProbeSets :Test(1) {
    my $self = shift;

    is_deeply($self->{fetched}->[0]->get_all_ProbeSets,
              $self->{expected}->{probe_sets},
              'get_all_ProbeSets() works');
}

sub get_array_chip_ids :Test(1) {
    my $self = shift;

    is_deeply($self->{fetched}->[0]->get_array_chip_ids,
              $self->{expected}->{array_chip_ids},
              'get_array_chip_ids() works');
}

sub get_ArrayChips :Test(1) {
    my $self = shift;

    is_deeply($self->{fetched}->[0]->get_ArrayChips,
              $self->{expected}->{array_chips},
              'get_ArrayChips() works');
}

sub get_design_ids :Test(1) {
    my $self = shift;

    is_deeply($self->{fetched}->[0]->get_design_ids,
              $self->{expected}->{design_ids},
              'get_design_ids() works');
}

1;