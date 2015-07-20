#!/usr/bin/env perl

use strict;
use warnings;
use Test::More;
use Test::Exception;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Funcgen::Channel;

# Module compiles
ok( 1, 'Bio::EnsEMBL::Funcgen::Channel compiles' );

# Test constructor without arguments
my $channel = Bio::EnsEMBL::Funcgen::Channel->new();
isa_ok ($channel, 'Bio::EnsEMBL::Funcgen::Channel', 'Channel');

# Test getters setters
ok( test_getter_setter( $channel, "sample_id", 'new_sample_id' ));
ok( test_getter_setter( $channel, "experimental_chip_id", 'new_ec_id' ));
ok( test_getter_setter( $channel, "type", 'new_type' ));
ok( test_getter_setter( $channel, "dye", 'new_dye' ));

# Test constructor with arguments
my $new_channel = Bio::EnsEMBL::Funcgen::Channel->new(
    -EXPERIMENTAL_CHIP_ID => 1234,
    -SAMPLE_ID            => 'sample id',
    -TYPE                 => 'type',
    -DYE                  => 'dye',
);

done_testing();
