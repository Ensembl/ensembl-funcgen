#!/usr/bin/env perl

use strict;
use warnings;
use Test::More;
use Test::Exception;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Funcgen::CoordSystem;
use Bio::EnsEMBL::Test::MultiTestDB;

# Module compiles
ok( 1, 'Bio::EnsEMBL::Funcgen::CoordSystem compiles' );

# Test constructor with mandatory arguments only
my $cosys = Bio::EnsEMBL::Funcgen::CoordSystem->new( -name => 'chromosome' );

isa_ok( $cosys, 'Bio::EnsEMBL::Funcgen::CoordSystem', 'CoordSystem' );

# Test name definition is present
throws_ok { Bio::EnsEMBL::Funcgen::CoordSystem->new }
qr/A name argument is required/, 'Check that name is supplied';

# Test constructor with optional arguments
$cosys = Bio::EnsEMBL::Funcgen::CoordSystem->new(
    -name    => 'chromosome',
    -version => 'NCBI33'
);

isa_ok( $cosys, 'Bio::EnsEMBL::Funcgen::CoordSystem', 'CoordSystem' );

# Test subroutines
$cosys->add_core_coord_system_info(
    -schema_build         => '53_34u',
    -core_coord_system_id => 2,
    -rank                 => 1,
);

is( $cosys->name(), 'chromosome', 'Retrieve coord system name' );

is( $cosys->version(), 'NCBI33', 'Retrieve coord system version' );

is( $cosys->get_latest_schema_build(), '53_34u', 'Retrieve schema build' );

is( $cosys->contains_schema_build('53_34u'),
    1, 'Contains schema build is true' );
is( $cosys->contains_schema_build('53_35u'),
    0, 'Contains schema build is false' );
throws_ok { $cosys->contains_schema_build() } qr/Must pass a schema_build/,
    'Check that schema build is supplied';

throws_ok { $cosys->is_top_level() }
qr/Not yet implmented, need to test against the core cache using dnadb\/schema_build/,
    'Check is_top_level throws exception';

# my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
# my $db   = $multi->get_DBAdaptor( "funcgen" );

# print $cosys->($db);

# 'SCHEMA_BUILD','TOP_LEVEL', 'SEQUENCE_LEVEL', 'DEFAULT', 'RANK', 'IS_STORED', 'CORE_COORD_SYSTEM_ID'

done_testing();
