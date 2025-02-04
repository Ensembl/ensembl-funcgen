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

package RegulatoryFeatureAdaptorTest;

use strict;
use warnings;

use Test::More;
use Test::Exception;

use Bio::EnsEMBL::Funcgen::RegulatoryFeature;
use Bio::EnsEMBL::Funcgen::RegulatoryActivity;
use parent qw(Bio::EnsEMBL::Funcgen::Test);

sub parameters :Test(setup) {
    my $self = shift;

    my $regulatory_build   = $self->_quick_fetch('RegulatoryBuild', 1);
    my $regulatory_feature = $self->_quick_fetch('RegulatoryFeature', 1);

    my $new_activity = Bio::EnsEMBL::Funcgen::RegulatoryActivity->new(
        '-activity'     => 'ACTIVE',
        '-epigenome_id' => 127
    );

    my $new_regulatory_feature =
        Bio::EnsEMBL::Funcgen::RegulatoryFeature->new_fast({
            slice                => $regulatory_feature->slice(),
            start                => 100_000,
            end                  => 110_000,
            strand               => 0,
            feature_type         => $regulatory_feature->get_FeatureType,
            _bound_lengths       => [ 50, 100 ],
            epigenome_count      => 20,
            stable_id            => 'ENSR00001000000',
            _analysis            => $regulatory_feature->get_Analysis,
            regulatory_build_id  => 1,
            _regulatory_activity => []
        });

    $new_regulatory_feature->add_regulatory_activity($new_activity);

    my $slice_adaptor = $self->{core_db}->get_adaptor('Slice');
    my $slice         = $slice_adaptor->fetch_by_seq_region_id(131545);

    my %parameters = (
        'stable_id'              => 'ENSR00001036347',
        'regulatory_build'       => $regulatory_build,
        'regulatory_feature'     => $regulatory_feature,
        'new_regulatory_feature' => $new_regulatory_feature,
        'slice'                  => $slice,
        'regulatory_build'       => $regulatory_build,
    );

    $self->{parameters} = \%parameters;
}

sub define_expected :Test(setup) {
    my $self = shift;

    my $regulatory_feature = $self->_quick_fetch('RegulatoryFeature', 1);
    my $valid_activities   =
        [ 'INACTIVE', 'REPRESSED', 'POISED', 'ACTIVE', 'NA' ];
    my $valid_activities_string = 'INACTIVE, REPRESSED, POISED, ACTIVE, NA';

    $self->{expected} = {
        'regulatory_feature'      => $regulatory_feature,
        'valid_activities'        => $valid_activities,
        'valid_activities_string' => $valid_activities_string,
        'rf_count'                => 1
    };

}

sub fetch_by_stable_id :Test(1) {
    my $self = shift;

    my $adaptor = $self->{funcgen_db}->get_adaptor('RegulatoryFeature');

    is_deeply(
        $adaptor->fetch_by_stable_id($self->{parameters}->{stable_id}),
        $self->{expected}->{regulatory_feature},
        'fetch_by_stable_id() works'
    );
}

sub fetch_by_stable_id_RegulatoryBuild :Test(1) {
    my $self = shift;

    my $adaptor = $self->{funcgen_db}->get_adaptor('RegulatoryFeature');

    my $stable_id        = $self->{parameters}->{stable_id};
    my $regulatory_build = $self->{parameters}->{regulatory_build};

    is_deeply(
        $adaptor->fetch_by_stable_id_RegulatoryBuild($stable_id,
                                                     $regulatory_build),
        $self->{expected}->{regulatory_feature},
        'fetch_by_stable_id_RegulatoryBuild() works'
    );
}

sub fetch_Iterator :Test(1) {
    my $self = shift;

    my $adaptor = $self->{funcgen_db}->get_adaptor('RegulatoryFeature');

    is_deeply(
        $adaptor->fetch_Iterator()->next(),
        $self->{expected}->{regulatory_feature},
        'fetch_Iterator() works'
    );
}

sub fetch_Iterator_by_RegulatoryBuild :Test(1) {
    my $self = shift;

    my $adaptor =
        $self->{funcgen_db}->get_adaptor('RegulatoryFeature');
    my $regulatory_build = $self->{parameters}->{regulatory_build};

    is_deeply(
        $adaptor->fetch_Iterator_by_RegulatoryBuild($regulatory_build)->next(),
        $self->{expected}->{regulatory_feature},
        'fetch_Iterator_by_RegulatoryBuild() works'
    );
}

sub store :Test(1) {
    my $self = shift;

    my $adaptor = $self->{funcgen_db}->get_adaptor('RegulatoryFeature');

    my $new_regulatory_feature = $self->{parameters}->{new_regulatory_feature};

    $self->{multi}->save('funcgen', 'regulatory_feature',
                         'regulatory_activity');

    $adaptor->store($new_regulatory_feature);

    my $newly_stored_rf =
        $adaptor->fetch_by_stable_id($new_regulatory_feature->stable_id);

    my $got = [
        $newly_stored_rf->seq_region_name(),
        $newly_stored_rf->seq_region_start(),
        $newly_stored_rf->seq_region_end(),
        $newly_stored_rf->bound_start_length(),
        $newly_stored_rf->bound_end_length(),
        $newly_stored_rf->seq_region_strand(),
        $newly_stored_rf->get_FeatureType(),
        $newly_stored_rf->stable_id(),
        $newly_stored_rf->epigenome_count(),
        $newly_stored_rf->regulatory_build_id(),
        $newly_stored_rf->regulatory_activity()->[0]->epigenome_id(),
        $newly_stored_rf->regulatory_activity()->[0]->activity(),
    ];

    my $expected = [
        $new_regulatory_feature->seq_region_name(),
        $new_regulatory_feature->seq_region_start(),
        $new_regulatory_feature->seq_region_end(),
        $new_regulatory_feature->bound_start_length(),
        $new_regulatory_feature->bound_end_length(),
        $new_regulatory_feature->seq_region_strand(),
        $new_regulatory_feature->get_FeatureType(),
        $new_regulatory_feature->stable_id(),
        $new_regulatory_feature->epigenome_count(),
        $new_regulatory_feature->regulatory_build_id(),
        $new_regulatory_feature->regulatory_activity()->[0]->epigenome_id(),
        $new_regulatory_feature->regulatory_activity()->[0]->activity(),
    ];

    is_deeply($got, $expected, 'store() method works');

    $self->{multi}->restore('funcgen', 'regulatory_feature',
                            'regulatory_activity');
}

sub valid_activities :Test(1) {
    my $self = shift;

    my $adaptor = $self->{funcgen_db}->get_adaptor('RegulatoryFeature');

    my @fetched_valid_activities = $adaptor->valid_activities;

    is_deeply(
        \@fetched_valid_activities,
        $self->{expected}->{valid_activities},
        'valid_activities() works'
    );
}

sub valid_activities_as_string :Test(1) {
    my $self = shift;

    my $adaptor = $self->{funcgen_db}->get_adaptor('RegulatoryFeature');

    is_deeply(
        $adaptor->valid_activities_as_string,
        $self->{expected}->{valid_activities_string},
        'valid_activities_as_string() works'
    );
}

sub fetch_all_by_Slice :Test(2) {
    my $self = shift;

    my $adaptor = $self->{funcgen_db}->get_adaptor('RegulatoryFeature');
    my $slice   = $self->{parameters}->{slice};

    is(
        scalar @{$adaptor->fetch_all_by_Slice($slice)},
        $self->{expected}->{rf_count},
        'fetch_all_by_Slice() returns the correct number of features'
    );

    is(
        $adaptor->fetch_all_by_Slice($slice)->[0]->stable_id(),
        $self->{expected}->{regulatory_feature}->stable_id(),
        'fetch_all_by_Slice() works'
    );
}

sub fetch_all_by_Slice_RegulatoryBuild :Test(2) {
    my $self = shift;

    my $adaptor = $self->{funcgen_db}->get_adaptor('RegulatoryFeature');

    my $slice            = $self->{parameters}->{slice};
    my $regulatory_build = $self->{parameters}->{regulatory_build};

    is(
        scalar @{$adaptor->fetch_all_by_Slice_RegulatoryBuild($slice,
                                                              $regulatory_build)},
        $self->{expected}->{rf_count},
        'fetch_all_by_Slice_RegulatoryBuild() returns the correct number of features'
    );

    is(
        $adaptor->fetch_all_by_Slice_RegulatoryBuild($slice, $regulatory_build)
                ->[0]->stable_id(),
        $self->{expected}->{regulatory_feature}->stable_id(),
        'fetch_all_by_Slice_RegulatoryBuild() works'
    );
}

sub fetch_all_by_Slice_FeatureSets :Test(1) {
    my $self = shift;

    my $adaptor = $self->{funcgen_db}->get_adaptor('RegulatoryFeature');

    dies_ok {
        $adaptor->fetch_all_by_Slice_FeatureSets();
    }
        'fetch_all_by_Slice_FeatureSets() dies as expected'
}

1;