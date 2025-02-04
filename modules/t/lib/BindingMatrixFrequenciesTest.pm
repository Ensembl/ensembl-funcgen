=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2025] EMBL-European Bioinformatics Institute

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

package BindingMatrixFrequenciesTest;

use strict;
use warnings;

use Test::More;
use Test::Exception;

use parent qw(Bio::EnsEMBL::Funcgen::Test);

sub parameters :Test(setup) {
    my $self = shift;

    my %mandatory_constructor_parameters = (
        -position   => 5,
        -nucleotide => 'G',
        -frequency  => 1000
    );

    $self->{mandatory_constructor_parameters} =
        \%mandatory_constructor_parameters;

    my $binding_matrix =
        $self->{funcgen_db}->get_adaptor('BindingMatrix')->fetch_by_dbID(1);

    my %optional_constructor_parameters = (
        -binding_matrix => $binding_matrix,
    );

    my %constructor_parameters = (%mandatory_constructor_parameters,
                                  %optional_constructor_parameters);

    $self->{constructor_parameters} = \%constructor_parameters;
}

sub define_expected :Test(setup) {
    my $self = shift;

    my $binding_matrix =
        $self->{funcgen_db}->get_adaptor('BindingMatrix')->fetch_by_dbID(134);

    $self->{expected} = {
        'binding_matrix' => $binding_matrix,
        'position'       => 1,
        'nucleotide'     => 'A',
        'frequency'      => 789
    };
}

sub dbIDs_to_fetch {return [8813];}

sub getters {
    return [ 'position', 'nucleotide', 'frequency' ];
}

sub get_BindingMatrix :Test(1) {
    my $self = shift;

    is_deeply($self->{fetched}->[0]->get_BindingMatrix(),
              $self->{expected}->{binding_matrix},
              'get_BindingMatrix works'
    );
}

1;
