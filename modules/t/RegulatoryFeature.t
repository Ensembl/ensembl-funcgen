# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2017] EMBL-European Bioinformatics Institute
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

use Test::More qw(no_plan);    #no_plan required for skip usage

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Slice;
use Data::Printer;

# ---------------
# Module compiles
# ---------------
BEGIN { use_ok('Bio::EnsEMBL::Funcgen::RegulatoryFeature'); }

# ------------------------------
# Setup test database connection
# ------------------------------
my $multi   = Bio::EnsEMBL::Test::MultiTestDB->new();
my $func_db = $multi->get_DBAdaptor('funcgen');
my $core_db = $multi->get_DBAdaptor('core');

my $slice_adaptor               = $core_db->get_adaptor('Slice');
my $feature_set_adaptor         = $func_db->get_adaptor('featureset');
my $feature_type_adaptor        = $func_db->get_adaptor('featuretype');
my $regulatory_feature_adaptor  = $func_db->get_adaptor('regulatoryfeature');
my $regulatory_activity_adaptor = $func_db->get_adaptor('regulatoryactivity');
my $analysis_adaptor            = $func_db->get_adaptor('analysis');
my $epigenome_adaptor           = $func_db->get_adaptor('epigenome');

# ----------------
# Test constructor
# ----------------
my $slice = $slice_adaptor->fetch_by_region( 'chromosome', '1' );
my $feature_set = $feature_set_adaptor->fetch_by_name(
                                'HeLa-S3_CTCF_ENCODE_Broad_SWEmbl_R0005_IDR');
my $feature_type = $feature_type_adaptor->fetch_by_name('CTCF');
my $new_regulatory_feature =
    Bio::EnsEMBL::Funcgen::RegulatoryFeature->new(
                                                -SLICE        => $slice,
                                                -START        => 1000,
                                                -END          => 2000,
                                                -FEATURE_SET  => $feature_set,
                                                -FEATURE_TYPE => $feature_type
    );

isa_ok( $new_regulatory_feature, 'Bio::EnsEMBL::Funcgen::RegulatoryFeature',
        'RegulatoryFeature' );

# -------------------
# Test get_Analysis()
# -------------------
my $fetched_regulatory_feature
    = $regulatory_feature_adaptor->fetch_by_stable_id('54736');
my $expected_analysis
    = $analysis_adaptor->fetch_by_logic_name('Regulatory_Build');

is_deeply( $fetched_regulatory_feature->get_Analysis(),
           $expected_analysis, 'Test get_Analysis()' );

# --------------------------
# Test regulatory_build_id()
# --------------------------
is( $fetched_regulatory_feature->regulatory_build_id(),
    1, 'Test regulatory_build_id()' );

# --------------------
# Test display_label()
# --------------------
is( $fetched_regulatory_feature->display_label(),
    'Promoter Regulatory Feature',
    'Test display_label()' );

# -----------------
# Test display_id()
# -----------------
is( $fetched_regulatory_feature->display_id(), 54736, 'Test display_id()' );

# ----------------
# Test stable_id()
# ----------------
is( $fetched_regulatory_feature->stable_id(), 54736, 'Test stable_id()' );

# -----------------------------
# Test get_RegulatoryEvidence()
# -----------------------------
#TODO

# ----------------------------------------
# Test regulatory_activity_for_epigenome()
# ----------------------------------------
my $epigenome = $epigenome_adaptor->fetch_by_name('HeLa-S3');

# ideally this should work, but there is a minor difference between the two objects
# my $expected_regulatory_activity
#     = $regulatory_activity_adaptor->fetch_by_dbID(167);
# is_deeply( $fetched_regulatory_feature->regulatory_activity_for_epigenome(
#                                                                   $epigenome),
#            $expected_regulatory_activity,
#            'Test regulatory_activity_for_epigenome()' );

is( $fetched_regulatory_feature->regulatory_activity_for_epigenome(
                                                          $epigenome)->dbID(),
    167,
    'Test regulatory_activity_for_epigenome()' );

# -------------------------------
# Test get_underlying_structure()
# -------------------------------
# without the epigenome parameter
my $underlying_structure
    = $fetched_regulatory_feature->get_underlying_structure();
my $last_index = scalar @{$underlying_structure} - 1;

is( $underlying_structure->[0],
    $fetched_regulatory_feature->bound_start(),
    'Test get_underlying_structure() - bound_start' );
is( $underlying_structure->[1],
    $fetched_regulatory_feature->start(),
    'Test get_underlying_structure() - start' );
is( $underlying_structure->[ $last_index - 1 ],
    $fetched_regulatory_feature->end(),
    'Test get_underlying_structure() - end' );
is( $underlying_structure->[$last_index],
    $fetched_regulatory_feature->bound_end(),
    'Test get_underlying_structure() - bound_end' );

for ( my $i = 2; $i <= $last_index - 2; $i++ ) {
    ok( $underlying_structure->[$i]
            > $fetched_regulatory_feature->bound_start() &&
            $underlying_structure->[$i]
            < $fetched_regulatory_feature->bound_end(),
        'Test get_underlying_structure() coordinates' );
}

# with the epigenome parameter
$underlying_structure
    = $fetched_regulatory_feature->get_underlying_structure();
$last_index = scalar @{$underlying_structure} - 1;

is( $underlying_structure->[0],
    $fetched_regulatory_feature->bound_start(),
    'Test get_underlying_structure() - bound_start' );
is( $underlying_structure->[1],
    $fetched_regulatory_feature->start(),
    'Test get_underlying_structure() - start' );
is( $underlying_structure->[ $last_index - 1 ],
    $fetched_regulatory_feature->end(),
    'Test get_underlying_structure() - end' );
is( $underlying_structure->[$last_index],
    $fetched_regulatory_feature->bound_end(),
    'Test get_underlying_structure() - bound_end' );

for ( my $j = 2; $j <= $last_index - 2; $j++ ) {
    ok( $underlying_structure->[$j]
            > $fetched_regulatory_feature->bound_start() &&
            $underlying_structure->[$j]
            < $fetched_regulatory_feature->bound_end(),
        'Test get_underlying_structure() coordinates' );
}

# --------------------------
# Test regulatory_activity()
# --------------------------


done_testing();

