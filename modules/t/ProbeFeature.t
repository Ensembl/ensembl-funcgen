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

# use diagnostics;
use autodie;

use Data::Dumper qw( Dumper );
use Test::More;
use Test::Exception;    # throws_ok
use Bio::EnsEMBL::Test::TestUtils qw( test_getter_setter debug );

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

# ---------------
# Module compiles
# ---------------
BEGIN { use_ok('Bio::EnsEMBL::Funcgen::ProbeFeature'); }

# ------------------------------
# Setup test database connection
# ------------------------------
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db    = $multi->get_DBAdaptor("funcgen");
my $pa    = $db->get_adaptor("probe");
my $sa    = $db->get_adaptor("slice");
my $pfa   = $db->get_adaptor("probe_feature");
my $psa   = $db->get_adaptor("probeset");

# ----------------
# Test constructor
# ----------------
# without arguments
my $probe_feature = Bio::EnsEMBL::Funcgen::ProbeFeature->new();

isa_ok( $probe_feature,
        'Bio::EnsEMBL::Funcgen::ProbeFeature',
        'ProbeFeature constructor return type (no arguments passed)' );

# with arguments
my $probe = $pa->fetch_by_array_probe_probeset_name( 'HumanWG_6_V2',
                                                     'ILMN_1677794' );
my $slice = $sa->fetch_by_region( 'chromosome', '1' );

$probe_feature =
    Bio::EnsEMBL::Funcgen::ProbeFeature->new( -PROBE         => $probe,
                                              -MISMATCHCOUNT => 0,
                                              -SLICE         => $slice,
                                              -START         => 10,
                                              -END           => 50,
                                              -STRAND        => 1,
                                              -ADAPTOR       => $pfa, );

isa_ok( $probe_feature,
        'Bio::EnsEMBL::Funcgen::ProbeFeature',
        'ProbeFeature constructor return type (arguments passed)' );

# fast version
#TODO: Test fast version of constructor

# ----------------------
# Test getters - setters
# ----------------------
my $probeset = $psa->fetch_all_by_name('214727_at')->[0];

ok( test_getter_setter( $probe_feature, 'probeset', $probeset ),
    'test_getter_setter ProbeFeature::probeset' );

#TODO: Test probeset_id()

ok( test_getter_setter( $probe_feature, 'mismatchcount', 5 ),
    'test_getter_setter ProbeFeature::mismatchcount' );

ok( test_getter_setter( $probe_feature, 'cigar_string', 'cigar_string' ),
    'test_getter_setter ProbeFeature::cigar_string' );

ok( test_getter_setter( $probe_feature, 'probe', $probe ),
    'test_getter_setter ProbeFeature::probe' );

my $not_a_probe_obj = $probeset;
throws_ok { $probe_feature->probe($not_a_probe_obj) }
qr /Probe must be a Bio::EnsEMBL::Funcgen::Probe object/,
    'Test probe() argument exception';

is( $probe_feature->probe_id(), 138152, 'Test probe_id()' );

# --------------------
# Test summary_as_hash
# --------------------
my $expected_summary = { array_probe_names => {
                                            'HumanHT-12_V3' => "ILMN_1677794",
                                            'HumanHT-12_V4' => "ILMN_1677794",
                                            'HumanRef-8_V3' => "ILMN_1677794",
                                            'HumanWG_6_V1'  => "0005670053",
                                            'HumanWG_6_V2'  => "ILMN_1677794",
                                            'HumanWG_6_V3'  => "ILMN_1677794"
                         },
                         end             => 50,
                         feature_type    => "array_probe",
                         probe_length    => 50,
                         seq_region_name => 1,
                         start           => 10,
                         strand          => 1 };

is_deeply( $probe_feature->summary_as_hash,
           $expected_summary, 'Test summary_as_hash()' );

done_testing();
