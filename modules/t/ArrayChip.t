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
BEGIN { use_ok('Bio::EnsEMBL::Funcgen::ArrayChip'); }

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

my $arraychip_adaptor = $db->get_ArrayChipAdaptor();

# ----------------
# Test constructor
# ----------------
my $array_chip = Bio::EnsEMBL::Funcgen::ArrayChip->new(
    -design_id => 'Test_design_id',
    -array_id  => 34,                  # HG-U133A
    -name      => 'Test_name',
    -adaptor   => $arraychip_adaptor
);                                     # Required for get_Array with -array_id

isa_ok(
    $array_chip,
    'Bio::EnsEMBL::Funcgen::ArrayChip',
    'ArrayChip constructor return type'
);

throws_ok {
    my $array_chip = Bio::EnsEMBL::Funcgen::ArrayChip->new(
        -design_id => 'Test_design_id',
        -array_id  => 34,                  # HG-U133A
        -adaptor   => $arraychip_adaptor
    );
}
qr/Must define a name(.*) and design_id(.*)/,
    'Test that a name and design_id is provided';

my $array = Bio::EnsEMBL::Funcgen::Array->new(
    -NAME        => 'Test',
    -FORMAT      => 'Tiled',
    -SIZE        => '1',
    -VENDOR      => 'VENDOR',
    -TYPE        => 'OLIGO',
    -DESCRIPTION => 'TESTING',
    -CLASS       => 'VENDOR_FORMAT',    #e.g. AFFY_UTR, ILLUMINA_WG
);

throws_ok {
    my $array_chip = Bio::EnsEMBL::Funcgen::ArrayChip->new(
    	-name      => 'Test_name',
        -design_id => 'Test_design_id',
        -array_id  => 34,                  # HG-U133A
        -array     => $array,
        -adaptor   => $arraychip_adaptor
    );
}
qr/Must provide either -array or -array_id but not both/,
    'Test that either an array or array_id is provided but not both';

$array = "not_a_Bio::EnsEMBL::Funcgen::Array";

throws_ok {
    my $array_chip = Bio::EnsEMBL::Funcgen::ArrayChip->new(
    	-name      => 'Test_name',
        -design_id => 'Test_design_id',
        -array     => $array,
        -adaptor   => $arraychip_adaptor
    );
}
qr/Array paramter must be a valid Bio::EnsEMBL::Funcgen::Array/,
    'Test that the array object provided is a Bio::EnsEMBL::Funcgen::Array';

# ----------------
# Test getters
# ----------------
is( $array_chip->name,      'Test_name',      'ArrayChip::name' );
is( $array_chip->design_id, 'Test_design_id', 'ArrayChip::design_id' );
is( $array_chip->array_id,  34,               'ArrayChip::array_id' );

# ----------------
# Test get_Array()
# ----------------
isa_ok(
    $array_chip->get_Array,
    'Bio::EnsEMBL::Funcgen::Array',
    'ArrayChip::get_Array return type'
);

done_testing();
