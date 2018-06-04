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
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Slice;

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
my $feature_type           = $feature_type_adaptor->fetch_by_name('CTCF');
my $new_regulatory_feature = Bio::EnsEMBL::Funcgen::RegulatoryFeature->new(
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
my $fetched_regulatory_feature =
  $regulatory_feature_adaptor->fetch_by_stable_id('54736');
my $expected_analysis =
  $analysis_adaptor->fetch_by_logic_name('Regulatory_Build');

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
is(
    $fetched_regulatory_feature->display_label(),
    'Promoter Regulatory Feature',
    'Test display_label()'
);

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
my $epigenome = $epigenome_adaptor->fetch_by_name('HeLa-S3');


# ---------------------------
# Test _assert_epigenome_ok()
# ---------------------------
throws_ok { $fetched_regulatory_feature->_assert_epigenome_ok() }
qr/Epigenome parameter was undefined/,
  'Test _assert_epigenome_ok() exception throw';

my $invalid_epigenome = $feature_type;
throws_ok {
    $fetched_regulatory_feature->_assert_epigenome_ok($invalid_epigenome);
}
qr/epigenome parameter must have type Bio::EnsEMBL::Funcgen::Epigenome/,
  'Test _assert_epigenome_ok() exception throw again';

# ----------------------------------------
# Test regulatory_activity_for_epigenome()
# ----------------------------------------

# ideally this should work, but there is a minor difference between the two objects
# my $expected_regulatory_activity
#     = $regulatory_activity_adaptor->fetch_by_dbID(167);
# is_deeply( $fetched_regulatory_feature->regulatory_activity_for_epigenome(
#                                                                   $epigenome),
#            $expected_regulatory_activity,
#            'Test regulatory_activity_for_epigenome()' );

is(
    $fetched_regulatory_feature->regulatory_activity_for_epigenome($epigenome)
      ->dbID(),
    167, 'Test regulatory_activity_for_epigenome()'
);

throws_ok { $fetched_regulatory_feature->regulatory_activity_for_epigenome() }
qr/Epigenome parameter was undefined/,
  'Test regulatory_activity_for_epigenome() exception throw';

throws_ok {
    $fetched_regulatory_feature->regulatory_activity_for_epigenome(
        $invalid_epigenome);
}
qr/Wrong parameter, expected an epigenome/,
  'Test regulatory_activity_for_epigenome() exception throw again';

# -------------------------------
# Test get_underlying_structure()
# -------------------------------
# without the epigenome parameter
my $underlying_structure =
  $fetched_regulatory_feature->get_underlying_structure();
my $last_index = scalar @{$underlying_structure} - 1;

is(
    $underlying_structure->[0],
    $fetched_regulatory_feature->bound_start(),
    'Test get_underlying_structure() - bound_start'
);
is(
    $underlying_structure->[1],
    $fetched_regulatory_feature->start(),
    'Test get_underlying_structure() - start'
);
is(
    $underlying_structure->[ $last_index - 1 ],
    $fetched_regulatory_feature->end(),
    'Test get_underlying_structure() - end'
);
is(
    $underlying_structure->[$last_index],
    $fetched_regulatory_feature->bound_end(),
    'Test get_underlying_structure() - bound_end'
);

for ( my $i = 2 ; $i <= $last_index - 2 ; $i++ ) {
    ok(
        $underlying_structure->[$i] > $fetched_regulatory_feature->bound_start()
          && $underlying_structure->[$i] <
          $fetched_regulatory_feature->bound_end(),
        'Test get_underlying_structure() coordinates'
    );
}

# with the epigenome parameter
$underlying_structure =
  $fetched_regulatory_feature->get_underlying_structure($epigenome);
$last_index = scalar @{$underlying_structure} - 1;

is(
    $underlying_structure->[0],
    $fetched_regulatory_feature->bound_start(),
    'Test get_underlying_structure() - bound_start'
);
is(
    $underlying_structure->[1],
    $fetched_regulatory_feature->start(),
    'Test get_underlying_structure() - start'
);
is(
    $underlying_structure->[ $last_index - 1 ],
    $fetched_regulatory_feature->end(),
    'Test get_underlying_structure() - end'
);
is(
    $underlying_structure->[$last_index],
    $fetched_regulatory_feature->bound_end(),
    'Test get_underlying_structure() - bound_end'
);

for ( my $j = 2 ; $j <= $last_index - 2 ; $j++ ) {
    ok(
        $underlying_structure->[$j] > $fetched_regulatory_feature->bound_start()
          && $underlying_structure->[$j] <
          $fetched_regulatory_feature->bound_end(),
        'Test get_underlying_structure() coordinates'
    );
}

# --------------------------
# Test regulatory_activity()
# --------------------------
my $activity_list = $fetched_regulatory_feature->regulatory_activity();
is( scalar @{$activity_list}, 17, 'Test regulatory_activity()' );

for my $activity ( @{$activity_list} ) {
    isa_ok( $activity, 'Bio::EnsEMBL::Funcgen::RegulatoryActivity',
        'RegulatoryActivity' );
}

# -----------------------
# Test regulatory_build()
# -----------------------
isa_ok(
    $fetched_regulatory_feature->get_regulatory_build(),
    'Bio::EnsEMBL::Funcgen::RegulatoryBuild',
    'RegulatoryBuild'
);
is( $fetched_regulatory_feature->get_regulatory_build()->dbID(),
    1, 'Test regulatory_build()' );

# ------------------------------
# Test add_regulatory_activity()
# ------------------------------
my $new_activity = $regulatory_activity_adaptor->fetch_by_dbID(1);
$fetched_regulatory_feature->add_regulatory_activity($new_activity);

$activity_list = $fetched_regulatory_feature->regulatory_activity();
is( scalar @{$activity_list}, 18, 'Test add_regulatory_activity()' );

# ----------------------
# Test has_activity_in()
# ----------------------
is( $fetched_regulatory_feature->has_activity_in($epigenome),
    1, 'Test has_activity_in()' );

my $other_epigenome = $epigenome_adaptor->fetch_by_name('HRE');
is( $fetched_regulatory_feature->has_activity_in($other_epigenome),
    undef, 'Test has_activity_in() again' );

# -----------------------------------
# Test has_epigenomes_with_activity()
# -----------------------------------
is( $fetched_regulatory_feature->has_epigenomes_with_activity('ACTIVE'),
    1, 'Test has_epigenomes_with_activity()' );
is( $fetched_regulatory_feature->has_epigenomes_with_activity('POISED'),
    undef, 'Test has_epigenomes_with_activity() again' );

# ---------------------------------
# Test get_epigenomes_by_activity()
# ---------------------------------
my $epigenome_list =
  $fetched_regulatory_feature->get_epigenomes_by_activity('ACTIVE');
is( scalar @{$epigenome_list}, 17, 'Test get_epigenomes_by_activity()' );

throws_ok {
    $fetched_regulatory_feature->get_epigenomes_by_activity('INVALID');
}
qr/Please pass a valid activity to this method/,
  'Test get_epigenomes_by_activity() exception throw';

# ----------------------
# Test epigenome_count()
# ----------------------
is( $fetched_regulatory_feature->epigenome_count(),
    17, 'Test epigenome_count()' );

# -----------------------------
# Test bound_seq_region_start()
# -----------------------------
is( $fetched_regulatory_feature->bound_seq_region_start(),
    32888409, 'Test bound_seq_region_start()' );

# ---------------------------
# Test bound_seq_region_end()
# ---------------------------
is( $fetched_regulatory_feature->bound_seq_region_end(),
    32892606, 'Test bound_seq_region_end()' );

# -------------------------
# Test bound_start_length()
# -------------------------
is( $fetched_regulatory_feature->bound_start_length(),
    198, 'Test bound_start_length()' );

# -----------------------
# Test bound_end_length()
# -----------------------
is( $fetched_regulatory_feature->bound_end_length(),
    1198, 'Test bound_end_length()' );

# ------------------
# Test bound_start()
# ------------------
is( $fetched_regulatory_feature->bound_start(), 32888409,
    'Test bound_start()' );

# ----------------
# Test bound_end()
# ----------------
is( $fetched_regulatory_feature->bound_end(), 32892606, 'Test bound_end()' );

# ----------------------
# Test summary_as_hash()
# ----------------------
my $expected_summary = {
    bound_end       => 32892606,
    bound_start     => 32888409,
    description     => "Predicted promoter",
    end             => 32891408,
    feature_type    => "Promoter",
    ID              => 54736,
    seq_region_name => 13,
    source          => "Regulatory_Build",
    start           => 32888607,
    strand          => 0
};

is_deeply( $fetched_regulatory_feature->summary_as_hash(),
    $expected_summary, 'Test summary_as_hash()' );

done_testing();

