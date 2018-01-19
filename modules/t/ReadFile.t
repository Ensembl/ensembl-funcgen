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
use autodie;
use feature qw(say);

# use Data::Dumper qw( Dumper );
use Test::More;
use Test::Exception;    # throws_ok
use Bio::EnsEMBL::Test::TestUtils qw( test_getter_setter debug );

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

# ---------------
# Module compiles
# ---------------
BEGIN { use_ok('Bio::EnsEMBL::Funcgen::ReadFile'); }

# ------------------------------
# Setup test database connection
# ------------------------------
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db    = $multi->get_DBAdaptor("funcgen");
my $epia   = $db->get_adaptor("epigenome");
my $ea    = $db->get_adaptor("experiment");
my $fta   = $db->get_adaptor("featuretype");
my $aa    = $db->get_adaptor("analysis");
my $issa  = $db->get_adaptor("ReadFile");

# ----------------
# Test constructor
# ----------------
my $epigenome    = $epia->fetch_by_name('HeLa-S3');
my $exp          = $ea->fetch_by_name('K562_WCE_ENCODE_UCHICAGO');
my $feature_type = $fta->fetch_by_name('CTCF');
my $analysis     = $aa->fetch_by_logic_name('FANTOM_v5');
my $is_control   = 1;

my $input_subset =
    Bio::EnsEMBL::Funcgen::ReadFile->new( -name => 'my_new_input_subset',
                                             -feature_type => $feature_type,
                                             -epigenome   => $epigenome,
                                             -experiment   => $exp,
                                             -analysis     => $analysis,
                                             -is_control   => $is_control,
                                             -replicate    => 1 );

isa_ok( $input_subset,
        'Bio::EnsEMBL::Funcgen::ReadFile',
        'ReadFile constructor return type' );

throws_ok {
    my $input_subset = Bio::EnsEMBL::Funcgen::ReadFile->new(
        -name         => 'my_new_input_subset',
        -feature_type => $feature_type,
        # -cell_type    => $cell_type,
        -experiment => $exp,
        -analysis   => $analysis,
        -is_control => $is_control, );
}
qr/Mandatory parameter -epigenome is not defined/,
    "Test constructor epigenome exception";

throws_ok {
    my $input_subset = Bio::EnsEMBL::Funcgen::ReadFile->new(
        -name         => 'my_new_input_subset',
        -feature_type => $feature_type,
        -epigenome    => $epigenome,
        # -experiment => $exp,
        -analysis   => $analysis,
        -is_control => $is_control, );
}
qr/Mandatory parameter -experiment is not defined/,
    "Test constructor experiment exception";

throws_ok {
    my $input_subset = Bio::EnsEMBL::Funcgen::ReadFile->new(
        -name         => 'my_new_input_subset',
        -feature_type => $feature_type,
        -epigenome    => $epigenome,
        -experiment   => $exp,
        -analysis     => $analysis,
        # -is_control => $is_control,
    );
}
qr/Must defined an -is_control parameter/,
    "Test constructor is_control exception";

# ------------
# Test getters
# ------------
# is( $input_subset->replicate(),  1, "Test ReadFile::replicate()" );
is( $input_subset->is_control(), 1, "Test ReadFile::is_control()" );

# ----------------------------------
# Test reset_relational_attributes()
# ----------------------------------
my $iss              = $issa->fetch_by_name('SRR037563');
my $expected_adaptor = $iss->{adaptor};

my $not_a_hashref = [ "not", "a", "hashref" ];

throws_ok {
    $input_subset->reset_relational_attributes($not_a_hashref);
}
qr/Must pass a HASHREF, not/,
    "Test reset_relational_attributes() parameter exception";

# test without resetting dbID and adaptor
$iss->reset_relational_attributes( {  -feature_type => $feature_type,
                                      -epigenome    => $epigenome,
                                      -experiment   => $exp,
                                      -analysis     => $analysis, },
                                   1 );

is( $iss->{feature_type}, $feature_type,
    "Test reset_relational_attributes() - feature_type" );
is( $iss->{epigenome}, $epigenome,
    "Test reset_relational_attributes() - cell_type" );
is( $iss->{experiment}, $exp,
    "Test reset_relational_attributes() - experiment" );
is( $iss->{analysis}, $analysis,
    "Test reset_relational_attributes() - analysis" );
is( $iss->{dbID}, 2366, "Test reset_relational_attributes() - dbID" );
is( $iss->{adaptor}, $expected_adaptor,
    "Test reset_relational_attributes() - adaptor" );

# now test after resetting dbID and adaptor
$iss->reset_relational_attributes( {  -feature_type => $feature_type,
                                      -epigenome    => $epigenome,
                                      -experiment   => $exp,
                                      -analysis     => $analysis, },
                                   0 );

is( $iss->{dbID}, undef, "Test reset_relational_attributes() - reset dbID" );
is( $iss->{adaptor}, undef,
    "Test reset_relational_attributes() - reset adaptor" );

# -----------------
# Test compare_to()
# -----------------
# my $new_iss = $issa->fetch_by_name('SRR037563');
# my $expected_comparison = { analysis     => "dbID mismatch:	71	-	61",
#                             cell_type    => "dbID mismatch:	1	-	5",
#                             experiment   => "dbID mismatch:	4	-	1041",
#                             feature_type => "dbID mismatch:	9	-	1",
#                             is_control   => "Return size mismatch:	1, 0",
#                             name         => "my_new_input_subset - SRR037563",
#                             replicate    => "1 - 3" };

# is_deeply( $input_subset->compare_to($new_iss),
#     $expected_comparison, "Test ReadFile::compare_to()" );

done_testing();
