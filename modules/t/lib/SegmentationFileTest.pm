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

package SegmentationFileTest;

use strict;
use warnings;

use Test::More;
use Test::Exception;

use parent qw(Bio::EnsEMBL::Funcgen::Test);

sub parameters :Test(setup) {
    my $self = shift;

    my %mandatory_constructor_parameters = ();

    $self->{mandatory_constructor_parameters} =
        \%mandatory_constructor_parameters;

    my $analysis         = $self->_quick_fetch('Analysis', 120);
    my $epigenome        = $self->_quick_fetch('Epigenome', 126);
    my $regulatory_build = $self->_quick_fetch('RegulatoryBuild', 1);
    my $segmentation     = $self->_quick_fetch('Segmentation', 1);

    my %optional_constructor_parameters = (
        dbID             => 123,
        name             => 'dummy_name',
        analysis         => $analysis,
        epigenome        => $epigenome,
        regulatory_build => $regulatory_build,
        file             => '/my/dummy/file',
        file_type        => 'FILE_TYPE',
        md5sum           => '1a2b3c4d5e',
        segmentation     => $segmentation
    );

    my %constructor_parameters = (%mandatory_constructor_parameters,
                                  %optional_constructor_parameters);

    $self->{constructor_parameters} = \%constructor_parameters;
}

sub define_expected :Test(setup) {
    my $self = shift;

    my $analysis         = $self->_quick_fetch('Analysis', 120);
    my $epigenome        = $self->_quick_fetch('Epigenome', 126);
    my $regulatory_build = $self->_quick_fetch('RegulatoryBuild', 1);
    my $segmentation_id  = 1;
    my $segmentation     =
        $self->_quick_fetch('Segmentation', $segmentation_id);

    my $adaptor = $self->{funcgen_db}->get_adaptor('Segmentation');

    $self->{expected} = {
        'dbID'              => 1,
        'name'              =>
        'Segmentation of neut. from V. in blueprint_no_ctcf',
        '_analysis'         => $analysis,
        '_epigenome'        => $epigenome,
        '_regulatory_build' => $regulatory_build,
        '_segmentation'     => $segmentation,
        'file'              => '/funcgen/segmentation_file/M1_CB.bb',
        'file_type'         => 'BIGBED',
        'adaptor'           => $adaptor,
        'md5sum'            => 'c142f681a6806fde296172c54ad03e5a',
        'segmentation_id'   => $segmentation_id,
    };
}

sub dbIDs_to_fetch {return [ 1 ];}

sub getters_setters {
    return [ 'dbID', 'name', '_analysis', '_epigenome', '_regulatory_build',
             '_segmentation', 'adaptor', 'md5sum', 'file', 'file_type',
             'segmentation_id' ];
}

sub get_Analysis :Test(1) {
    my $self = shift;

    is_deeply($self->{fetched}->[0]->get_Analysis,
              $self->{expected}->{_analysis},
              'get_Analysis() works'
    );
}

sub get_Epigenome :Test(1) {
    my $self = shift;

    is_deeply($self->{fetched}->[0]->get_Epigenome,
              $self->{expected}->{_epigenome},
              'get_Epigenome() works'
    );
}

sub get_RegulatoryBuild :Test(1) {
    my $self = shift;

    is_deeply($self->{fetched}->[0]->get_RegulatoryBuild,
              $self->{expected}->{_regulatory_build},
              'get_RegulatoryBuild() works'
    );
}

sub get_Segmentation :Test(1) {
    my $self = shift;

    is_deeply($self->{fetched}->[0]->get_Segmentation,
              $self->{expected}->{_segmentation},
              'get_Segmentation() works'
    );
}

1;

