# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#!/usr/bin/env perl

use strict;
use warnings;
use Test::More;
use Test::Exception;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Funcgen::Channel;

# Module compiles
use_ok('Bio::EnsEMBL::Funcgen::Channel');

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
