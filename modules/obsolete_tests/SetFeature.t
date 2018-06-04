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
BEGIN { use_ok('Bio::EnsEMBL::Funcgen::SetFeature'); }

# ------------------------------
# Setup test database connection
# ------------------------------
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db    = $multi->get_DBAdaptor("funcgen");
my $fsa   = $db->get_adaptor("featureset");
my $fta   = $db->get_adaptor("featuretype");
my $aa    = $db->get_adaptor("analysis");
my $epia   = $db->get_adaptor("epigenome");

# ----------------
# Test constructor
# ----------------
my $feature_set
    = $fsa->fetch_by_name('A549_CTCF_ENCODE_Broad_SWEmbl_R0005_IDR');
my $ft       = $fta->fetch_by_name('CTCF');
my $analysis = $aa->fetch_by_logic_name('SWEmbl_R0005_IDR');

my $set_feature =
    Bio::EnsEMBL::Funcgen::SetFeature->new( -SET           => $feature_set,
                                            -FEATURE_TYPE  => $ft,
                                            -ANALYSIS      => $analysis,
                                            -DISPLAY_LABEL => 'test_label', );

isa_ok( $set_feature,
        'Bio::EnsEMBL::Funcgen::SetFeature',
        'SetFeature constructor return type' );

throws_ok {
    my $set_feature =
        Bio::EnsEMBL::Funcgen::SetFeature->new(
                                  -SET => 'not a FeatureSet/ResultSet object',
                                  -FEATURE_TYPE  => $ft,
                                  -ANALYSIS      => $analysis,
                                  -DISPLAY_LABEL => 'test_label', );
}
qr/Must pass valid Bio::EnsEMBL::Funcgen::FeatureSet or ResultSet object/,
    "Constructor's exception for FeatureSet/ResultSet parameter";

throws_ok {
    my $set_feature =
        Bio::EnsEMBL::Funcgen::SetFeature->new(
                                  -SET          => $feature_set,
                                  -FEATURE_TYPE => 'not a FeatureType object',
                                  -ANALYSIS     => $analysis,
                                  -DISPLAY_LABEL => 'test_label', );
}
qr/feature_type param must be a valid Bio::EnsEMBL::Funcgen::FeatureType/,
    "Constructor's exception for FeatureType parameter";

# --------------
# Test getters()
# --------------
is( $set_feature->feature_set(),
    $feature_set, 'Test SetFeature::set_feature() getter' );
is( $set_feature->set(), $feature_set, 'Test SetFeature::set() getter' );

my $epigenome = $epia->fetch_by_name('A549');
is_deeply( $set_feature->epigenome(),
           $epigenome, 'Test SetFeature::epigenome() getter' );

is( $set_feature->feature_type(),
    $ft, 'Test SetFeature::feature_type() getter' );

my $new_set_feature = Bio::EnsEMBL::Funcgen::SetFeature->new(
    -SET => $feature_set,
    # -FEATURE_TYPE  => $ft,
    -ANALYSIS      => $analysis,
    -DISPLAY_LABEL => 'test_label', );
is_deeply( $new_set_feature->feature_type(),
    $ft,
    'Test SetFeature::feature_type() getter without predefined FeatureType' );

is_deeply( $set_feature->analysis(),
    $analysis, 'Test SetFeature::analysis() getter' );
is( $set_feature->display_label(),
    'test_label', 'Test SetFeature::display_label() getter' );

done_testing();
