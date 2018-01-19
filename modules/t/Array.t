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

use Bio::EnsEMBL::Funcgen::ArrayChip;
use Bio::EnsEMBL::Funcgen::Array;

use Test::More;
use Test::Exception;    # throws_ok
use Bio::EnsEMBL::Test::TestUtils qw( test_getter_setter debug );
use Bio::EnsEMBL::Test::MultiTestDB;

# ---------------
# Module compiles
# ---------------
BEGIN { use_ok('Bio::EnsEMBL::Funcgen::Array'); }

# ------------------------------
# Setup test database connection
# ------------------------------
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db    = $multi->get_DBAdaptor("funcgen");

isa_ok(
    $db,
    'Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor',
    'Test database instatiated'
);

my $array_adaptor = $db->get_ArrayAdaptor();

# ----------------
# Test constructor
# ----------------
my $array = Bio::EnsEMBL::Funcgen::Array->new(
    -NAME        => 'Test',
    -FORMAT      => 'Tiled',
    -SIZE        => '1',
    -VENDOR      => 'VENDOR',
    -TYPE        => 'OLIGO',
    -DESCRIPTION => 'TESTING',
    -CLASS       => 'VENDOR_FORMAT',    #e.g. AFFY_UTR, ILLUMINA_WG
);

isa_ok(
    $array,
    'Bio::EnsEMBL::Funcgen::Array',
    'Array constructor return type'
);

throws_ok {
    my $array = Bio::EnsEMBL::Funcgen::Array->new(
        -NAME        => 'Test',
        -FORMAT      => 'Tiled',
        -SIZE        => '1',
        -TYPE        => 'OLIGO',
        -DESCRIPTION => 'TESTING',
        -CLASS       => 'VENDOR_FORMAT',    #e.g. AFFY_UTR, ILLUMINA_WG
    );
}
qr/Must provide a vendor parameter/, 'Test that a vendor is provided';

throws_ok {
    my $array = Bio::EnsEMBL::Funcgen::Array->new(
        -FORMAT      => 'Tiled',
        -SIZE        => '1',
        -VENDOR      => 'VENDOR',
        -TYPE        => 'OLIGO',
        -DESCRIPTION => 'TESTING',
        -CLASS       => 'VENDOR_FORMAT',    #e.g. AFFY_UTR, ILLUMINA_WG
    );
}
qr/Must provide a name parameter/, 'Test that a name is provided';

# ----------------------
# Test getters - setters
# ----------------------
ok( test_getter_setter( $array, 'name', 'Test2' ),
    'test_getter_setter Array::name' );
ok( test_getter_setter( $array, 'format', 'TILED2' ),
    'test_getter_setter Array::format' );
ok( test_getter_setter( $array, 'vendor', 'ILLUMINA2' ),
    'test_getter_setter Array::vendor' );
ok( test_getter_setter( $array, 'class', 'ILLUMINA_WG2' ),
    'test_getter_setter Array::class' );
ok( test_getter_setter( $array, 'description', 'TESTING2' ),
    'test_getter_setter Array::description' );
ok( test_getter_setter( $array, 'type', 'OLIGO2' ),
    'test_getter_setter Array::type' );

# ---------------------
# Test get_all_Probes()
# ---------------------
$array = $array_adaptor->fetch_by_name_vendor( 'SurePrint_G3_GE_8x60k',
    'AGILENT' );

my $probe_adaptor   = $db->get_ProbeAdaptor;
my $expected_probes = $probe_adaptor->fetch_all_by_Array($array);

is_deeply( $array->get_all_Probes(),
    $expected_probes, 'Test get_all_Probes() subroutine' );

## --------------------------
## Test get_all_Probe_dbIDs()
## --------------------------
#$array = $array_adaptor->fetch_by_name_vendor( 'HG-U133A', 'AFFY' );
#
#my @expected_probe_dbIDs = (
#    '655066',  '655071',  '655076',  '655080',  '655085',  '655090',
#    '655095',  '655100',  '655104',  '655109',  '655114',  '1967802',
#    '2481435', '5724311', '5724315', '5724318', '5724322', '5724325',
#    '5724328', '5724332', '5724335', '5724338', '5724341', '5724345'
#);
#
#is( @{ $array->get_all_Probe_dbIDs() },
#    @expected_probe_dbIDs, 'Test get_all_Probe_dbIDs() subroutine' );
#
#throws_ok {
#    my $array = Bio::EnsEMBL::Funcgen::Array->new(
#        -NAME        => 'Test',
#        -FORMAT      => 'Tiled',
#        -SIZE        => '1',
#        -VENDOR      => 'VENDOR',
#        -TYPE        => 'OLIGO',
#        -DESCRIPTION => 'TESTING',
#        -CLASS       => 'VENDOR_FORMAT',    #e.g. AFFY_UTR, ILLUMINA_WG
#    );
#    $array->get_all_Probe_dbIDs();
#}
#qr/Must have set an adaptor to get_all_Probe_dbIDs/,
#    'Test get_all_Probe_dbIDs(), adaptor is provided';

$array = $array_adaptor->fetch_by_name_vendor( 'HG-U133A', 'AFFY' );

# ------------------------
# Test get_all_ProbeSets()
# ------------------------
my $probeset_adaptor   = $db->get_ProbeSetAdaptor();
my $expected_probesets = $probeset_adaptor->fetch_all_by_Array($array);

is_deeply( $array->get_all_ProbeSets(),
    $expected_probesets, 'Test get_all_ProbeSets() subroutine' );

# -------------------------
# Test get_array_chip_ids()
# -------------------------
is_deeply( $array->get_array_chip_ids(),
    [74], 'Test get_array_chip_ids() subroutine' );

# ---------------------
# Test get_design_ids()
# ---------------------
is_deeply( $array->get_design_ids(),
    ['HG-U133A'], 'Test get_design_ids() subroutine' );

# ------------------
# Test probe_count()
# ------------------
is( $array->probe_count(), 24, 'Test probe_count() subroutine' );

# ---------------------------------
# Test get_ArrayChip_by_design_id()
# ---------------------------------
my $expected_array_chip
    = $db->get_ArrayChipAdaptor->fetch_by_array_design_ids( 34, 'HG-U133A' );
is_deeply( $array->get_ArrayChip_by_design_id('HG-U133A'),
    $expected_array_chip, 'Test get_ArrayChip_by_design_id() subroutine' );

done_testing();