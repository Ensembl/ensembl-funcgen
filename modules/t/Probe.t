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

use Test::More;
use Test::Exception;    # throws_ok
use Test::Warn;
use Bio::EnsEMBL::Test::TestUtils qw( test_getter_setter debug );

use Bio::EnsEMBL::Test::MultiTestDB;
use Data::Dumper;
use feature qw(say);

# ---------------
# Module compiles
# ---------------
BEGIN { use_ok('Bio::EnsEMBL::Funcgen::Probe'); }

# ------------------------------
# Setup test database connection
# ------------------------------
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db    = $multi->get_DBAdaptor("funcgen");
my $pa    = $db->get_adaptor("probe");
my $aa    = $db->get_adaptor("array");
my $psa   = $db->get_adaptor("probeset");

# ----------------------------------------------------
# Test constructor, using -sequence and -probe_seq_id
# ----------------------------------------------------
my $array = $aa->fetch_by_name_vendor( 'HuEx-1_0-st-v2', 'AFFY' );
my $probe_set = $psa->fetch_all_by_name('2715818');

my $probe = Bio::EnsEMBL::Funcgen::Probe->new(
    -NAME          => 'probe_with_sequence',
    -ARRAY         => $array,
    -PROBE_SET     => $probe_set,
    -CLASS         => "EXPERIMENTAL",
    -SEQUENCE      => 'AGAGTATAGGTGTGATGATG',
);
isa_ok(
    $probe,
    'Bio::EnsEMBL::Funcgen::Probe',
    'Probe constructor return type [sequence]'
);

$probe = Bio::EnsEMBL::Funcgen::Probe->new(
    -NAME          => 'probe_with_probe_seq_id' ,
    -ARRAY         => $array,
    -PROBE_SET     => $probe_set,
    -CLASS         => "EXPERIMENTAL",
    -PROBE_SEQ_ID  => 1,
    );
isa_ok(
    $probe,
    'Bio::EnsEMBL::Funcgen::Probe',
    'Probe constructor return type [probe_seq_id]'
);
#---------------------------------------------------
# Test constructor throws
#---------------------------------------------------
throws_ok {
my $probe = Bio::EnsEMBL::Funcgen::Probe->new(
    -NAME          => 'probe_with_no_probe_seq_id_and_no_sequence' ,
    -ARRAY         => $array,
    -PROBE_SET     => $probe_set,
    -CLASS         => "EXPERIMENTAL",
    );
}
qr /Define either -SEQUENCE or PROBE_SEQ_ID/,
    'Sequence information missing';

throws_ok {
my $probe = Bio::EnsEMBL::Funcgen::Probe->new(
    -NAME          => 'probe_with_no_probe_seq_id_and_no_sequence' ,
    -ARRAY         => $array,
    -PROBE_SET     => $probe_set,
    -CLASS         => "EXPERIMENTAL",
    -SEQUENCE      => 'AGAGTATAGGTGTGATGATG',
    -PROBE_SEQ_ID  => 1,
    );
}
qr /Either define -SEQUENCE or PROBE_SEQ_ID, not both/,
    'Sequence information duplicate';

throws_ok {
my $probe = Bio::EnsEMBL::Funcgen::Probe->new(
    -NAME          => 'probe_with_no_array' ,
    -PROBE_SET     => $probe_set,
    -CLASS         => "EXPERIMENTAL",
    -SEQUENCE      => 'AGAGTATAGGTGTGATGATG'
    );
}
qr /Bio::EnsEMBL::Funcgen::Array required/,
    'Probe must be attached to an Array';

throws_ok {
my $probe = Bio::EnsEMBL::Funcgen::Probe->new(
    -NAME          => 'probe_with_object_that_is_not_an_array' ,
    -PROBE_SET     => $probe_set,
    -CLASS         => "EXPERIMENTAL",
    -SEQUENCE      => 'AGAGTATAGGTGTGATGATG',
    -ARRAY         => 'FAKE_ARRAY',
    );
}
qr /Not a Bio::EnsEMBL::Funcgen::Array object/,
    'Not an Array object';


# -------------------------------
# Test add_array_chip_probename()
# -------------------------------
$probe->add_array_chip_probename( 'probe_with_probe_seq_id', $array );

is_deeply(
    $probe->{probenames},
    {   'HuEx-1_0-st-v2' =>[
            'probe_with_probe_seq_id'
            ]

    },
    'Test add_array_chip_probename()'
);
# Not available in e89
# is_deeply(
#     $probe->{arrays},
#     { '68' => $array },
#     'Test add_array_chip_probename() again'
# );
#
my $not_an_array;
throws_ok {
    $probe->add_array_chip_probename( 'new_probe', $not_an_array );
}
qr/You must pass a valid Bio::EnsEMBL::Funcgen::Array/,
    'Test add_array_chip_probename() exception';
#print Dumper($probe);
## ----------------------------
## Test get_all_ProbeFeatures()
## ----------------------------
my $new_probe = $pa->fetch_by_array_probe_probeset_name( 'HumanWG_6_V2', 'ILMN_1677794' );

my $expected_features
    = $db->get_ProbeFeatureAdaptor->fetch_all_by_Probe($new_probe);

is_deeply( $new_probe->get_all_ProbeFeatures(),
    $expected_features, 'Test get_all_ProbeFeatures()' );
#
# ---------------------
# Test get_all_Arrays()
# ---------------------
my $expected_arrays = [  $new_probe->{array} ];
is_deeply( $new_probe->get_all_Arrays, $expected_arrays,
    'Test get_all_Arrays()' );

# -----------------------
# Test get_names_Arrays()
# -----------------------
my $expected_name_array_pairs = {$new_probe->{array}->name,$new_probe->{array}};
#say Dumper($new_probe->get_names_Arrays);
is_deeply( $new_probe->get_names_Arrays,
    $expected_name_array_pairs, 'Test get_names_Arrays()' );
#
## -------------------------
## Test get_all_probenames()
## -------------------------
my $expected_probenames = [ 'ILMN_1677794' ];
#
#say Dumper($new_probe);
is( @{ $new_probe->get_all_probenames },
    @{$expected_probenames}, 'Test get_all_probenames()' );
#
## --------------------
## Test get_probename()
## --------------------
is( $new_probe->get_probename('HumanWG_6_V2'),
    'ILMN_1677794', 'Test get_probename()' );
#
##TODO: More tests needed
#
# -----------------------------
# Test get_all_complete_names()
# -----------------------------
my $expected_complete_names = [ 'HumanWG_6_V2:ILMN_1677794' ];
is( @{ $new_probe->get_all_complete_names },
    @{$expected_complete_names},
    'Test get_all_complete_names()'
);
#
## ------------------------
## Test get_complete_name()
## ------------------------
throws_ok { $new_probe->get_complete_name() }
qr/Must provide and array name argument to retreive the complete name/,
    'Test get_complete_name() no argument exception';

throws_ok { $new_probe->get_complete_name('Invalid_array_name') }
qr/Probe linked to /,
    'Test get_complete_name() unknown array name exception';

#TODO: More tests needed, not just for exceptions

# ----------------------
# Test getters - setters
# ----------------------
my $probeset = $psa->fetch_all_by_Array($array)->[0];

ok( test_getter_setter( $new_probe, 'probeset', $probeset ),
    'test_getter_setter Probe::probeset' );

ok( test_getter_setter( $new_probe, 'class', 'EXPERIMENTAL' ),
    'test_getter_setter Probe::class' );

ok( test_getter_setter( $new_probe, 'length', 25 ),
    'test_getter_setter Probe::length' );

ok( test_getter_setter( $new_probe, 'description', 'Excellent probe' ),
    'test_getter_setter Probe::description' );

# --------------------
# Test feature_count()
# --------------------
# is( $new_probe->feature_count(), 1, 'Test feature_count()' );

done_testing();
