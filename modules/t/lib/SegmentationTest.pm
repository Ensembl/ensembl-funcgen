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

=cut

package SegmentationTest;

use strict;
use warnings;

use Test::More;
use Test::Exception;

use parent qw(Bio::EnsEMBL::Funcgen::Test);

sub define_expected :Test(setup) {
    my $self = shift;

    my $segmentation_adaptor = $self->{funcgen_db}->get_adaptor('Segmentation');

    $self->{expected} = {
        'dbID'                => 1,
        'name'                => 'encode_ctcf',
        'class'               => 'ctcf',
        'superclass'          => 'encode',
        'regulatory_build_id' => undef,
        'adaptor'             => $segmentation_adaptor,
        'db'                  => $segmentation_adaptor
    };
}

sub dbIDs_to_fetch {return [ 1 ];}

sub getters_setters {
    return [ 'dbID', 'name', 'class', 'superclass', 'regulatory_build_id',
             'adaptor', 'db' ]
}

1;