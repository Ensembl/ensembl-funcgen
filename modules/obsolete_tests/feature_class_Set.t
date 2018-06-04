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
use Test::Exception;    # throws_ok
use Test::Warn;
use Bio::EnsEMBL::Test::TestUtils qw( test_getter_setter debug );
use Data::Dumper;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Funcgen::ResultSet;

# ---------------
# Module compiles
# ---------------
BEGIN { use_ok('Bio::EnsEMBL::Funcgen::feature_class_Set'); }

# ------------------------------
# Setup test database connection
# ------------------------------
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db    = $multi->get_DBAdaptor("funcgen");
my $fta   = $db->get_adaptor("featuretype");
my $aa    = $db->get_adaptor("analysis");
my $rsa   = $db->get_adaptor("resultset");

# ------------------------------
# Test _validate_feature_class()
# ------------------------------
my $feature_type = $fta->fetch_by_name("H3K27ac");
my $analysis     = $aa->fetch_by_logic_name("ChiP-Seq");

my @constructor_args = (
    -name          => 'placeholder',
    -table_name    => 'input_subset',
    -analysis      => $analysis,
    -feature_class => 'result',
    -feature_type  => $feature_type,
    -adaptor       => $rsa,
);

my $result_set = Bio::EnsEMBL::Funcgen::ResultSet->new(@constructor_args);

my $feature_class
    = Bio::EnsEMBL::Funcgen::feature_class_Set::_validate_feature_class(
    $result_set, \@constructor_args );

is( $feature_class, 'result', 'Test _validate_feature_class() return value' );

# --------------------
# Test feature_class()
# --------------------
is( Bio::EnsEMBL::Funcgen::feature_class_Set::feature_class($result_set),
    'result', 'Test feature_class()' );

# -------------------------
# Test feature_class_name()
# -------------------------
is( Bio::EnsEMBL::Funcgen::feature_class_Set::feature_class_name($result_set),
    'ResultFeature', 'Test feature_class_name()'
);

done_testing();

