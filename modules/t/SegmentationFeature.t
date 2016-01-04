# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.



use strict;
use warnings;
use diagnostics;
use autodie;
use feature qw(say);

# use Data::Dumper qw( Dumper );
use Test::More;
use Test::Exception;    # throws_ok
use Bio::EnsEMBL::Test::TestUtils qw( test_getter_setter debug );

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

# ---------------
# Module compiles
# ---------------
BEGIN { use_ok('Bio::EnsEMBL::Funcgen::SegmentationFeature'); }

# ------------------------------
# Setup test database connection
# ------------------------------
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db    = $multi->get_DBAdaptor("funcgen");
my $fsa   = $db->get_adaptor("featureset");
my $fta   = $db->get_adaptor("featuretype");
my $sfa   = $db->get_adaptor("segmentationfeature");
my $sa    = $db->get_adaptor("slice");

# ----------------
# Test constructor
# ----------------
my $feature_set
    = $fsa->fetch_by_name('A549_CTCF_ENCODE_Broad_SWEmbl_R0005_IDR');
my $feature_type = $fta->fetch_by_name('H3K4me3');
my $slice = $sa->fetch_by_region( 'chromosome', 'X' );

my $segmentation_feature =
    Bio::EnsEMBL::Funcgen::SegmentationFeature->new(
                                               -FEATURE_SET  => $feature_set,
                                               -FEATURE_TYPE => $feature_type,
                                               -SCORE        => 10,
                                               -DISPLAY_LABEL => 'test_label',
                                               -START         => 1,
                                               -END           => 100,
                                               -SLICE         => $slice, );

isa_ok( $segmentation_feature,
        'Bio::EnsEMBL::Funcgen::SegmentationFeature',
        'SegmentationFeature constructor return type' );

throws_ok {
    my $segmentation_feature
        = Bio::EnsEMBL::Funcgen::SegmentationFeature->new(
        -FEATURE_SET => $feature_set,
        # -FEATURE_TYPE => $feature_type
        );
}
qr/You must pass a valid FeatureType/,
    "Test constructor's exception for mandatory feature_type parameter";

# ------------
# Test getters
# ------------
is( $segmentation_feature->score(), 10, 'Test SegmentationFeature::score()' );
is( $segmentation_feature->display_label(),
    'test_label', 'Test SegmentationFeature::display_label()' );

my $new_sf = Bio::EnsEMBL::Funcgen::SegmentationFeature->new(
    -FEATURE_SET  => $feature_set,
    -FEATURE_TYPE => $feature_type,
    -SCORE        => 10,
    -START        => 1,
    -END          => 100,
    -SLICE        => $slice,
    # -DISPLAY_LABEL => 'test_label',
    -ADAPTOR => $sfa, );

is( $new_sf->display_label(),
    'H3K4me3 - A549',
    'Test SegmentationFeature::display_label() without predefined label' );

my $expected_summary = { segmentation_feature_type => 'H3K4me3',
                         cell_type                 => 'A549',
                         start                     => 1,
                         end                       => 100,
                         seq_region_name           => 'X' };

is_deeply( $segmentation_feature->summary_as_hash(),
           $expected_summary,
           'Test SegmentationFeature::expected_summary()' );

done_testing();
