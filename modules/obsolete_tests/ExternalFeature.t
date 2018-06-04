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
use diagnostics;

use Data::Dumper qw( Dumper );
use Test::More;
use Test::Exception;    # throws_ok
use Bio::EnsEMBL::Test::TestUtils qw( test_getter_setter debug );

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis;

# ---------------
# Module compiles
# ---------------
BEGIN { use_ok('Bio::EnsEMBL::Funcgen::ExternalFeature'); }

# ------------------------------
# Setup test database connection
# ------------------------------
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db    = $multi->get_DBAdaptor("funcgen");
my $fta   = $db->get_adaptor("featuretype");
my $fsa   = $db->get_adaptor("featureset");
# my $ana   = $db->get_adaptor("analysis");
my $sla   = $db->get_adaptor("slice");
my $efa   = $db->get_adaptor("externalfeature");

# ----------------
# Test constructor
# ----------------
my $feature_set
    = $fsa->fetch_by_name('A549_CTCF_ENCODE_Broad_SWEmbl_R0005_IDR');
my $feature_type = $fta->fetch_by_name('CTCF');
# my $analysis     = $ana->fetch_by_logic_name('SWEmbl_R0005_IDR');
my $slice = $sla->fetch_by_name('chromosome:NCBI34:X:1000000:2000000:1');

my $external_feature = Bio::EnsEMBL::Funcgen::ExternalFeature->new(
    -FEATURE_SET => $feature_set,
    -SLICE       => $slice,
    -START       => 10,
    -END         => 100,
    -STRAND      => 1,
    -ADAPTOR     => $efa,
);

isa_ok(
    $external_feature,
    'Bio::EnsEMBL::Funcgen::ExternalFeature',
    'ExternalFeature constructor return type'
);

# --------------------
# Test display_label()
# --------------------
is($external_feature->display_label(),'CTCF - A549', 'Test display_label()');

done_testing();
