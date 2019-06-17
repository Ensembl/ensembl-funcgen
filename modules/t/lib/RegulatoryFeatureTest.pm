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

package RegulatoryFeatureTest;

use strict;
use warnings;

use Test::More;
use Test::Exception;

use Bio::EnsEMBL::Funcgen::RegulatoryActivity;
use parent qw(Bio::EnsEMBL::Funcgen::Test);

sub parameters :Test(setup) {
    my $self = shift;

    my %mandatory_constructor_parameters = ();

    $self->{mandatory_constructor_parameters} =
        \%mandatory_constructor_parameters;

    my $analysis = $self->_quick_fetch('Analysis', 16);

    my %optional_constructor_parameters = (
        -stable_id       => 'ENSR00000000000',
        -projected       => 1,
        -epigenome_count => 10,
        -analysis        => $analysis,
    );

    my %constructor_parameters = (%mandatory_constructor_parameters,
                                  %optional_constructor_parameters);

    $self->{constructor_parameters} = \%constructor_parameters;

    my $epigenome             = $self->_quick_fetch('Epigenome', 146);
    my $no_activity_epigenome = $self->_quick_fetch('Epigenome', 126);
    my $regulatory_activity   = $self->_quick_fetch('RegulatoryActivity', 64);

    my $additional_activity = Bio::EnsEMBL::Funcgen::RegulatoryActivity->new();

    my %parameters = (
        'epigenome'             => $epigenome,
        'no_activity_epigenome' => $no_activity_epigenome,
        'regulatory_activity'   => $regulatory_activity,
        'additional_activity'   => $additional_activity,
        'invalid_activity'      => 'INVALID'
    );

    $self->{parameters} = \%parameters;
}

sub define_expected :Test(setup) {
    my $self = shift;

    my $analysis       = $self->_quick_fetch('Analysis', 16);
    my $feature_type   = $self->_quick_fetch('FeatureType', 179363);
    my $display_label  = 'Promoter Flanking Region Regulatory Feature';
    my $motif_features =
        $self->_quick_fetch_all('MotifFeature',
                                [ 12383983, 12952012, 13512233 ]);
    my $verified_motif_feature = $self->_quick_fetch('MotifFeature', 9607738);
    my $regulatory_activity    =
        $self->_quick_fetch_all('RegulatoryActivity', [ 1, 18, 64 ]);
    my $regulatory_build = $self->_quick_fetch('RegulatoryBuild', 1);
    my $epigenome        = $self->_quick_fetch('Epigenome', 146);
    my $single_activity  = $self->_quick_fetch('RegulatoryActivity', 64);
    $single_activity->{epigenome} = $epigenome;

    $self->{expected} = {
        'analysis'               => $analysis,
        'feature_type'           => $feature_type,
        'regulatory_build_id'    => 1,
        'display_label'          => $display_label,
        'stable_id'              => 'ENSR00001036347',
        'single_activity'        => $single_activity,
        'motif_features'         => $motif_features,
        'verified_motif_feature' => $verified_motif_feature,
        'regulatory_activity'    => $regulatory_activity,
        'regulatory_build'       => $regulatory_build,
        'epigenome'              => $epigenome,
        'epigenome_count'        => 8,
        'bound_seq_region_start' => 113379701,
        'bound_seq_region_end'   => 113380400,
        'bound_start_length'     => 100,
        'bound_end_length'       => 200,
        'bound_start'            => 113379701,
        'bound_end'              => 113380400
    };

    my %summary = (
        'id'              => 'ENSR00001036347',
        'source'          => 'Regulatory_Build',
        'bound_start'     => 113379701,
        'bound_end'       => 113380400,
        'start'           => 113379801,
        'end'             => 113380200,
        'strand'          => 0,
        'seq_region_name' => 2,
        'description'     => 'Predicted promoter flanking region',
        'feature_type'    => 'Promoter Flanking Region'
    );

    $self->{expected}->{summary} = \%summary;
}

sub dbIDs_to_fetch {return [ 1, 64576 ];}

sub getters {
    return [ 'analysis', 'feature_type', 'regulatory_build_id',
             'display_label', 'stable_id', 'regulatory_activity',
             'epigenome_count', 'bound_seq_region_start',
             'bound_seq_region_end', 'bound_start_length', 'bound_end_length',
             'bound_start', 'bound_end' ];
}

sub get_Analysis :Test(1) {
    my $self = shift;

    is_deeply($self->{fetched}->[0]->get_Analysis(),
              $self->{expected}->{analysis},
              'get_Analysis() works'
    );
}

sub get_FeatureType :Test(1) {
    my $self = shift;

    is_deeply($self->{fetched}->[0]->get_FeatureType(),
              $self->{expected}->{feature_type},
              'get_FeatureType() works'
    );
}

sub display_id :Test(1) {
    my $self = shift;

    is($self->{fetched}->[0]->display_id(),
       $self->{expected}->{stable_id},
       'display_id() works'
    );
}

sub get_RegulatoryEvidence :Test(no_plan) {
    #TODO write tests
}

sub regulatory_activity_for_epigenome :Test(3) {
    my $self = shift;

    my $regulatory_feature = $self->{fetched}->[0];
    my $epigenome          = $self->{parameters}->{epigenome};
    my $not_an_epigenome   = $self->{parameters}->{regulatory_activity};

    is_deeply(
        $regulatory_feature->regulatory_activity_for_epigenome($epigenome),
        $self->{expected}->{single_activity},
        'regulatory_activity_for_epigenome() works'
    );
    # use Data::Printer;
    # p $regulatory_feature->regulatory_activity_for_epigenome($epigenome);
    # p $self->{expected}->{single_activity};
    # p $regulatory_feature->regulatory_activity;
    throws_ok {
        $regulatory_feature->regulatory_activity_for_epigenome();
    }
        qr/Epigenome parameter was undefined/,
        "... and exception is thrown when no parameter is passed";

    throws_ok {
        $regulatory_feature
            ->regulatory_activity_for_epigenome($not_an_epigenome);
    }
        qr/Wrong parameter, expected an epigenome, but got a/,
        "... and exception is thrown when the passed parameter is not an Epigenome object";
}

sub get_all_MotifFeatures :Test(1) {
    my $self = shift;

    is_deeply($self->{fetched}->[0]->get_all_MotifFeatures(),
              $self->{expected}->{motif_features},
              'get_all_MotifFeatures() works'
    );
}

sub get_all_experimentally_verified_MotifFeatures :Test(1) {
    my $self = shift;

    is_deeply(
        $self->{fetched}->[1]->get_all_experimentally_verified_MotifFeatures()
             ->[0],
        $self->{expected}->{verified_motif_feature},
        'get_all_experimentally_verified_MotifFeatures() works'
    );
}

sub get_regulatory_build :Test(1) {
    my $self = shift;

    is_deeply(
        $self->{fetched}->[0]->get_regulatory_build(),
        $self->{expected}->{regulatory_build},
        'get_regulatory_build() works'
    );
}

sub add_regulatory_activity :Test(1) {
    my $self = shift;

    $self->{fetched}->[0]->add_regulatory_activity($self->{parameters}
                                                        ->{additional_activity});

    is_deeply(
        pop $self->{fetched}->[0]->regulatory_activity(),
        $self->parameters->{additional_activity},
        'add_regulatory_activity() works'
    );
}

sub has_activity_in :Test(2) {
    my $self = shift;

    my $epigenome = $self->{parameters}->{epigenome};

    ok($self->{fetched}->[0]->has_activity_in($epigenome),
       'has_activity_in() works'
    );

    my $no_activity_epigenome = $self->{parameters}->{no_activity_epigenome};

    ok(!$self->{fetched}->[0]->has_activity_in($no_activity_epigenome),
       'has_activity_in() works again'
    );
}

sub has_epigenomes_with_activity :Test(2) {
    my $self = shift;

    my $activity = $self->{parameters}->{regulatory_activity}->activity;

    ok($self->{fetched}->[0]->has_epigenomes_with_activity($activity),
       'has_epigenomes_with_activity() works'
    );

    my $no_epigenome_activity =
        $self->{parameters}->{additional_activity}->activity;

    ok(!$self->{fetched}->[0]
             ->has_epigenomes_with_activity($no_epigenome_activity),
       'has_epigenomes_with_activity() works again'
    );
}

sub get_epigenomes_by_activity :Test(2) {
    my $self = shift;

    my $activity = $self->{parameters}->{regulatory_activity}->activity;

    is_deeply(
        $self->{fetched}->[0]->get_epigenomes_by_activity($activity),
        [ $self->{expected}->{epigenome} ],
        'get_epigenomes_by_activity() works'
    );

    my $invalid_activity = $self->{parameters}->{invalid_activity};

    throws_ok {
        $self->{fetched}->[0]->get_epigenomes_by_activity($invalid_activity)
    }
        qr/Please pass a valid activity to this method/,
        "... and exception is thrown when the passed activity is invalid";
}

1;
