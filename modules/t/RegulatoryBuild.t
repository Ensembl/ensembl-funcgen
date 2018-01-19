# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

use Test::More;
use Test::Exception;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;

# ---------------
# Module compiles
# ---------------
BEGIN { use_ok('Bio::EnsEMBL::Funcgen::RegulatoryBuild'); }

# ------------------------------
# Setup test database connection
# ------------------------------
my $multi   = Bio::EnsEMBL::Test::MultiTestDB->new();
my $func_db = $multi->get_DBAdaptor('funcgen');

my $analysis_adaptor           = $func_db->get_adaptor('analysis');
my $feature_type_adaptor       = $func_db->get_adaptor('featuretype');
my $regulatory_build_adaptor   = $func_db->get_adaptor('regulatorybuild');
my $regulatory_feature_adaptor = $func_db->get_adaptor('regulatoryfeature');

# ----------------
# Test constructor
# ----------------
my $new_regulatory_build = Bio::EnsEMBL::Funcgen::RegulatoryBuild->new(
    -name
    ,                       => 'The Ensembl Regulatory Build',
    -version                => '14',
    -initial_release_date   => '2016-6',
    -last_annotation_update => '2016-6',
    -feature_type_id        => 18,
    -analysis_id            => 15,
    -is_current             => 1
);

isa_ok(
    $new_regulatory_build,
    'Bio::EnsEMBL::Funcgen::RegulatoryBuild',
    'Regulatory Build'
);

# ----------------------
# Test getters - setters
# ----------------------
ok( test_getter_setter( $new_regulatory_build, 'name', 'renamed_RB' ),
    'test_getter_setter RegulatoryBuild::name()' );

ok(
    test_getter_setter( $new_regulatory_build, 'version', '15' ),
    'test_getter_setter RegulatoryBuild::version()'
);

ok(
    test_getter_setter(
        $new_regulatory_build, 'initial_release_date', '2016-7'
    ),
    'test_getter_setter RegulatoryBuild::initial_release_date()'
);

ok(
    test_getter_setter(
        $new_regulatory_build, 'last_annotation_update', '2016-7'
    ),
    'test_getter_setter RegulatoryBuild::last_annotation_update()'
);

ok(
    test_getter_setter( $new_regulatory_build, 'feature_type_id', '10' ),
    'test_getter_setter RegulatoryBuild::feature_type_id()'
);

ok(
    test_getter_setter( $new_regulatory_build, 'analysis_id', '10' ),
    'test_getter_setter RegulatoryBuild::analysis_id()'
);

ok(
    test_getter_setter( $new_regulatory_build, 'is_current', '0' ),
    'test_getter_setter RegulatoryBuild::is_current()'
);

ok(
    test_getter_setter(
        $new_regulatory_build, 'sample_regulatory_feature_id', '123'
    ),
    'test_getter_setter RegulatoryBuild::sample_regulatory_feature_id()'
);

# -------------------
# Test set_Analysis()
# -------------------
my $fetched_regulatory_build =
  $regulatory_build_adaptor->fetch_by_name('Regulatory features');

my $analysis = $analysis_adaptor->fetch_by_logic_name('Regulatory_Build');
$fetched_regulatory_build->set_Analysis($analysis);

is( $fetched_regulatory_build->analysis_id(), 16, 'Test set_Analysis()' );

# throws_ok { $fetched_regulatory_build->set_Analysis(); }
# qr /Analysis was not defined/,
#   'Test set_Analysis() exception throw';

# ---------------------
# Test fetch_Analysis()
# ---------------------
is_deeply( $fetched_regulatory_build->fetch_Analysis(),
    $analysis, 'Test fetch_Analysis()' );

# ----------------------
# Test set_FeatureType()
# ----------------------
my $feature_type = $feature_type_adaptor->fetch_by_name('RegulatoryFeature');
$fetched_regulatory_build->set_FeatureType($feature_type);

is( $fetched_regulatory_build->feature_type_id(), 19,
    'Test set_FeatureType()' );

# ------------------------
# Test fetch_FeatureType()
# ------------------------
is_deeply( $fetched_regulatory_build->fetch_FeatureType(),
    $feature_type, 'Test fetch_FeatureType()' );

# -------------------------
# Test get_all_Epigenomes()
# -------------------------

my $used_epigenomes = $fetched_regulatory_build->get_all_Epigenomes();

is( scalar @$used_epigenomes, 17, 'Test get_all_Epigenomes()' );

for my $used_epigenome ( @{$used_epigenomes} ) {
    isa_ok( $used_epigenome, 'Bio::EnsEMBL::Funcgen::Epigenome', 'Epigenome' );
}

# -----------------------------
# Test _get_all_epigenome_ids()
# -----------------------------
my $epigenome_ids = $fetched_regulatory_build->_get_all_epigenome_ids();
is( scalar @$epigenome_ids, 17, 'Test _get_all_epigenome_ids()' );

# -------------------------------------
# Test fetch_sample_RegulatoryFeature()
# -------------------------------------
# my $expected_sample_regulatory_feature =
#   $regulatory_feature_adaptor->fetch_by_stable_id('1639593');

# is_deeply(
#     $fetched_regulatory_build->fetch_sample_RegulatoryFeature(),
#     $expected_sample_regulatory_feature,
#     'Test fetch_sample_RegulatoryFeature()'
# );

done_testing();
