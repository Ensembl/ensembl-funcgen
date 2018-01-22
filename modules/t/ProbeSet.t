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
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Funcgen::ProbeSet;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Funcgen::DBSQL::ProbeSetAdaptor;

# ---------------
# Module compiles
# ---------------
BEGIN { use_ok('Bio::EnsEMBL::Funcgen::ProbeSet'); }

# ----------------------------------------------
# Test constructor with mandatory arguments only
# ----------------------------------------------
my $ps = Bio::EnsEMBL::Funcgen::ProbeSet->new();

isa_ok( $ps, 'Bio::EnsEMBL::Funcgen::ProbeSet', 'ProbeSet' );

# ----------------------------------------
# Test constructor with optional arguments
# ----------------------------------------
$ps = Bio::EnsEMBL::Funcgen::ProbeSet->new(
    -NAME   => 'ProbeSet-1',
    -SIZE   => 1,
    -FAMILY => "ENCODE_REGIONS",
);

isa_ok( $ps, 'Bio::EnsEMBL::Funcgen::ProbeSet', 'ProbeSet' );

# --------------------
# Test getters setters
# --------------------
ok( test_getter_setter( $ps, 'name',   'rc_AI071570_at' ), 'Test name() subroutine' );
ok( test_getter_setter( $ps, 'family', 'family' ), 'Test family() subroutine' );
ok( test_getter_setter( $ps, 'size',   '10' ), 'Test size() subroutine');

# ------------------------------
# Setup test database connection
# ------------------------------
my $multi      = Bio::EnsEMBL::Test::MultiTestDB->new();
my $func_db    = $multi->get_DBAdaptor("funcgen");
my $ps_adaptor = Bio::EnsEMBL::Funcgen::DBSQL::ProbeSetAdaptor->new($func_db);

# ------------------------------
# Test get_all_Arrays subroutine
# ------------------------------
$ps = $ps_adaptor->fetch_by_array_probe_set_name( 'PrimeView', '11740416_at' );
my $arrays = $ps->get_all_Arrays();

foreach my $array ( @{$arrays} ) {
    isa_ok( $array, 'Bio::EnsEMBL::Funcgen::Array' );
}

# ------------------------------
# Test get_all_Probes subroutine
# ------------------------------
my $probes = $ps->get_all_Probes();

foreach my $probe ( @{$probes} ) {
    isa_ok( $probe, 'Bio::EnsEMBL::Funcgen::Probe' );
}

done_testing();
