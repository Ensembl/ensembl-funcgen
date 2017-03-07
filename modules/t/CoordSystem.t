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
use Test::Exception;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Funcgen::CoordSystem;
use Bio::EnsEMBL::Funcgen::DBSQL::CoordSystemAdaptor;
use Bio::EnsEMBL::Funcgen::CellType;

# ---------------
# Module compiles
# ---------------
BEGIN { use_ok('Bio::EnsEMBL::Funcgen::CoordSystem'); }

# ----------------------------------------------
# Test constructor with mandatory arguments only
# ----------------------------------------------
my $cosys = Bio::EnsEMBL::Funcgen::CoordSystem->new( -name => 'chromosome' );

isa_ok( $cosys, 'Bio::EnsEMBL::Funcgen::CoordSystem', 'CoordSystem' );

# Test name definition is present
throws_ok { Bio::EnsEMBL::Funcgen::CoordSystem->new }
qr/A name argument is required/, 'Check that name is supplied';

# ----------------------------------------
# Test constructor with optional arguments
# ----------------------------------------
$cosys = Bio::EnsEMBL::Funcgen::CoordSystem->new(
    -name    => 'chromosome',
    -version => 'NCBI33'
);

isa_ok( $cosys, 'Bio::EnsEMBL::Funcgen::CoordSystem', 'CoordSystem' );

# ---------------------------------------------------------------------
# Test name, version, get_latest_schema_build, is_top_level subroutines
# ---------------------------------------------------------------------
$cosys->add_core_coord_system_info(
    -schema_build         => '53_34u',
    -core_coord_system_id => 2,
    -rank                 => 1,
);

is( $cosys->name(), 'chromosome', 'Retrieve coord system name' );

is( $cosys->version(), 'NCBI33', 'Retrieve coord system version' );

is( $cosys->get_latest_schema_build(), '53_34u', 'Retrieve schema build' );

throws_ok { $cosys->is_top_level() }
qr/Not yet implemented, need to test against the core cache using dnadb\/schema_build/,
    'Check is_top_level throws exception';

# -------------------------------------
# Test contains_schema_build subroutine
# -------------------------------------
is( $cosys->contains_schema_build('53_34u'),
    1, 'Contains schema build is true' );
is( $cosys->contains_schema_build('53_35u'),
    0, 'Contains schema build is false' );
throws_ok { $cosys->contains_schema_build() } qr/Must pass a schema_build/,
    'Check that schema build is supplied';

# ------------------------------------------
# Test add_core_coord_system_info subroutine
# ------------------------------------------
throws_ok {
    $cosys->add_core_coord_system_info(
        -core_coord_system_id => 2,
        -rank                 => 1,
    );
}
qr/Must provide a schema_build/, 'Test schema_build exception';

throws_ok {
    $cosys->add_core_coord_system_info(
        -schema_build => '53_34u',
        -rank         => 1,
    );
}
qr/Must provide a core_coord_system_id/,
    'Test core_coord_system_id exception';

throws_ok {
    $cosys->add_core_coord_system_info(
        -schema_build         => '53_34u',
        -core_coord_system_id => 2,
        -rank                 => 1,
        -top_level            => 1,
    );
}
qr/RANK argument must be 0 if TOP_LEVEL is 1/,
    'Test rank - top_level exception';

throws_ok {
    $cosys->add_core_coord_system_info(
        -schema_build         => '53_34u',
        -core_coord_system_id => 2,
        -rank                 => -1,
    );
}
qr/The RANK argument must be a positive integer/,
    'Test rank is positive integer exception';

throws_ok {
    $cosys->add_core_coord_system_info(
        -schema_build         => '53_34u',
        -core_coord_system_id => 2,
        -top_level            => 1,
    );
}
qr/The NAME argument must be "toplevel" if TOP_LEVEL is 1/,
    'Test name - top_level exception';

$cosys = Bio::EnsEMBL::Funcgen::CoordSystem->new(
    -name    => 'toplevel',
    -version => 'NCBI33'
);

throws_ok {
    $cosys->add_core_coord_system_info(
        -schema_build         => '53_34u',
        -core_coord_system_id => 2,
        -top_level            => 1,
        -sequence_level       => 2,
    );
}
qr/SEQUENCE_LEVEL argument must be 0 if TOP_LEVEL is 1/,
    'Test sequence_level - top_level exception';

throws_ok {
    $cosys->add_core_coord_system_info(
        -schema_build         => '53_34u',
        -core_coord_system_id => 2,
    );
}
qr/RANK argument must be non-zero if not toplevel CoordSystem/,
    'Test rank - top_level exception';

throws_ok {
    $cosys->add_core_coord_system_info(
        -schema_build         => '53_34u',
        -core_coord_system_id => 2,
        -rank                 => 2,
    );
}
qr/Cannot name coord system 'toplevel' unless TOP_LEVEL is 1/,
    'Test name - top_level exception';

# ------------------------------
# Setup test database connection
# ------------------------------
my $multi   = Bio::EnsEMBL::Test::MultiTestDB->new();

my $func_db = $multi->get_DBAdaptor("funcgen");
my $cs_func_adaptor
    = Bio::EnsEMBL::Funcgen::DBSQL::CoordSystemAdaptor->new($func_db);

# my $core_db = $multi->get_DBAdaptor("core");
# my $cs_core_adaptor
#     = Bio::EnsEMBL::Funcgen::DBSQL::CoordSystemAdaptor->new($core_db);

# ----------------------
# Test equals subroutine
# ----------------------
throws_ok { $cosys->equals(); }
qr/Argument must be a Bio::EnsEMBL::Funcgen::CoordSystem/,
    'Test argument is supplied for equals subroutine';

my $not_a_cosys = Bio::EnsEMBL::Funcgen::CellType->new( -name => "U2OS" );

throws_ok { $cosys->equals($not_a_cosys); }
qr/Argument must be a Bio::EnsEMBL::Funcgen::CoordSystem/,
    'Test the supplied argument for equals subroutine is a CoordSystem';

$cosys = Bio::EnsEMBL::Funcgen::CoordSystem->new(
    -name    => 'chromosome',
    -version => 'GRCh37',
    -adaptor => $cs_func_adaptor,
);

$cosys->add_core_coord_system_info(
    -schema_build         => '79_37',
    -core_coord_system_id => 2,
    -rank                 => 1
);

my $equal_cosys = $cosys;

is( $cosys->equals($equal_cosys), 1, 'Test two equal CoordSystem objects' );

my $not_equal_cosys = Bio::EnsEMBL::Funcgen::CoordSystem->new(
    -name    => 'contig',
    -version => '',
);

is( $cosys->equals($not_equal_cosys),
    0, 'Test two unequal CoordSystem objects' );

#test lines 436-446

# my $new_cosys =  Bio::EnsEMBL::Funcgen::CoordSystem->new(
#     -name    => 'chromosome',
#     -version => 'GRCh37',
#     -adaptor => $cs_core_adaptor,
# );

# $cosys->add_core_coord_system_info(
#     -schema_build         => '',
#     -core_coord_system_id => 2,
#     -rank                 => 1
# );
# # throws_ok{$cosys->equals($new_cosys);}
# # qr//,
# # '';
# $cosys->equals($new_cosys);

# ----------------------------------------------------
# Test rank, is_sequence_level, is_default subroutines
# ----------------------------------------------------
is( $cosys->rank($func_db), 1, 'Test rank subroutine' );
is( $cosys->is_sequence_level($func_db),
    0, 'Test is_sequence_level is false' );
is( $cosys->is_default($func_db), 0, 'Test is_default is false' );

# ------------------------------------------
# Test get_coord_system_attribute subroutine
# ------------------------------------------
$cosys->add_core_coord_system_info(
    -schema_build         => '79_37',
    -core_coord_system_id => 2,
    -rank                 => 1,
    -sequence_level       => 1,
    -default              => 1,
);

is( $cosys->is_sequence_level($func_db), 1,
    'Test is_sequence_level is true' );
is( $cosys->is_default($func_db), 1, 'Test is_default is true' );

throws_ok { $cosys->get_coord_system_attribute('is_stored'); }
qr/You must pass a dnadb to access the CoordSystem attribute/,
    'Test dnadb argument is supplied for get_coord_system_attribute subroutine';

my $not_a_DBAdaptor = Bio::EnsEMBL::Funcgen::CellType->new( -name => "U2OS" );

throws_ok {
    $cosys->get_coord_system_attribute( 'is_stored', $not_a_DBAdaptor );
}
qr/You must pass a dnadb to access the CoordSystem attribute/,
    'Test the supplied argument for get_coord_system_attribute subroutine is a Bio::EnsEMBL::DBSQL::DBAdaptor';

$cosys = Bio::EnsEMBL::Funcgen::CoordSystem->new(
    -name    => 'chromosome',
    -version => 'GRCh37',
    -adaptor => $cs_func_adaptor,
);
$cosys->add_core_coord_system_info(
    -schema_build         => 'non_existant',
    -core_coord_system_id => 2,
    -rank                 => 1
);

throws_ok { $cosys->get_coord_system_attribute( 'is_stored', $func_db ); }
qr/CoordSystem does not contain the schema_build/,
    'Test if CoordSystem contains the schema build for get_coord_system_attribute subroutine';

done_testing();
