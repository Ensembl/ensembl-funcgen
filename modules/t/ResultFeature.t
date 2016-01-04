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

use Test::More;
use Test::Exception;    # throws_ok
use Bio::EnsEMBL::Test::TestUtils qw( test_getter_setter debug );

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

# ---------------
# Module compiles
# ---------------
BEGIN { use_ok('Bio::EnsEMBL::Funcgen::ResultFeature'); }

# ------------------------------
# Setup test database connection
# ------------------------------
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db    = $multi->get_DBAdaptor("funcgen");
my $pa    = $db->get_adaptor("probe");
my $sa = $db->get_adaptor("slice");

# ----------------
# Test constructor
# ----------------
my $start  = 1;
my $end    = 100;
my $strand = -1;
my $score  = 5;
my $probe  = $pa->fetch_by_array_probe_probeset_name( 'HumanWG_6_V2',
                                                     'ILMN_1677794' );
my $result_set_id = 982;
my $window_size   = 6;

my $rf =
    Bio::EnsEMBL::Funcgen::ResultFeature->new_fast( $start, $end, $strand,
                               $score, $probe, $result_set_id, $window_size );
isa_ok( $rf,
        'Bio::EnsEMBL::Funcgen::ResultFeature',
        'ResultFeature constructor return type' );

# ------------
# Test getters
# ------------
my $slice = $sa->fetch_by_region( 'chromosome', 'X' );

is( $rf->start(), $start,
    'Test Bio::EnsEMBL::Funcgen::ResultFeature::start()' );
is( $rf->end(), $end, 'Test Bio::EnsEMBL::Funcgen::ResultFeature::end()' );
is( $rf->score(), $score,
    'Test Bio::EnsEMBL::Funcgen::ResultFeature::score()' );
is( $rf->strand(), $strand,
    'Test Bio::EnsEMBL::Funcgen::ResultFeature::strand()' );
is( $rf->probe(), $probe,
    'Test Bio::EnsEMBL::Funcgen::ResultFeature::probe()' );
is( $rf->result_set_id(), $result_set_id,
    'Test Bio::EnsEMBL::Funcgen::ResultFeature::result_set_id()' );
is( $rf->window_size(), $window_size,
    'Test Bio::EnsEMBL::Funcgen::ResultFeature::window_size()' );
is( $rf->slice($slice), $slice,
    'Test Bio::EnsEMBL::Funcgen::ResultFeature::slice()' );
is( $rf->length, 100, 'Test Bio::EnsEMBL::Funcgen::ResultFeature::length()' );

# -----------
# Test move()
# -----------
throws_ok {
    $rf->move();
}
qr/start and end arguments are required/,
    "Test Bio::EnsEMBL::Funcgen::ResultFeature::move(), presence of start and end arguments";

throws_ok {
    $rf->move( 50, 10 );
}
qr/start must be less than or equal to end/,
    "Test Bio::EnsEMBL::Funcgen::ResultFeature::move(), start is less or equal to end";

throws_ok {
    $rf->move( 5, 10, 3 );
}
qr/strand must be 0, -1 or 1/,
    "Test Bio::EnsEMBL::Funcgen::ResultFeature::move(), strand must be 0, -1 or 1";

$rf->move( 50, 200, 1 );

is( $rf->start(), 50,
    'Test Bio::EnsEMBL::Funcgen::ResultFeature::move(), start' );
is( $rf->end(), 200,
    'Test Bio::EnsEMBL::Funcgen::ResultFeature::move(), end' );
is( $rf->strand(), 1,
    'Test Bio::EnsEMBL::Funcgen::ResultFeature::move(), strand' );

# --------------------
# Test feature_Slice()
# --------------------
my $new_slice = $rf->feature_Slice();

isa_ok( $new_slice, 'Bio::EnsEMBL::Slice',
        'Bio::EnsEMBL::Funcgen::ResultFeature::feature_Slice return type' );
is( $new_slice->start(), 50,
    'Test Bio::EnsEMBL::Funcgen::ResultFeature::feature_Slice, start' );
is( $new_slice->end(), 200,
    'Test Bio::EnsEMBL::Funcgen::ResultFeature::feature_Slice, end' );
is( $new_slice->strand(), 1,
    'Test Bio::EnsEMBL::Funcgen::ResultFeature::feature_Slice, strand' );

done_testing();
