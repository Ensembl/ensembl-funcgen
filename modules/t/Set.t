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
BEGIN { use_ok('Bio::EnsEMBL::Funcgen::Set'); }

# ------------------------------
# Setup test database connection
# ------------------------------
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db    = $multi->get_DBAdaptor("funcgen");
my $fta   = $db->get_adaptor("featuretype");
my $aa    = $db->get_adaptor("analysis");
my $epia   = $db->get_adaptor("epigenome");
my $ea    = $db->get_adaptor("experiment");
my $issa  = $db->get_adaptor("inputsubset");

# ----------------
# Test constructor
# ----------------
my $feature_type = $fta->fetch_by_name('CTCF');
my $analysis     = $aa->fetch_by_logic_name('FANTOM_v5');
my $epigenome    = $epia->fetch_by_name('HeLa-S3');
my $exp          = $ea->fetch_by_name('K562_WCE_ENCODE_UCHICAGO');

my $set =
    Bio::EnsEMBL::Funcgen::Set->new( -name         => 'set_name',
                                     -feature_type => $feature_type,
                                     -analysis     => $analysis,
                                     -epigenome    => $epigenome,
                                     -dbID         => 1000,
                                     -adaptor      => $issa,
                                     -experiment   => $exp, );

isa_ok( $set, 'Bio::EnsEMBL::Funcgen::Set', 'Set constructor return type' );

throws_ok {
    my $set = Bio::EnsEMBL::Funcgen::Set->new(
        # -name         => 'set_name',
        -feature_type => $feature_type,
        -analysis     => $analysis,
        -epigenome    => $epigenome,
        -dbID         => 1000,
        -adaptor      => $issa,
        -experiment   => $exp, );
}
qr/Need to specify a name/, "Test constructor's name exception";

throws_ok {
    my $set = Bio::EnsEMBL::Funcgen::Set->new(
        -name => 'set_name',
        # -feature_type => $feature_type,
        -analysis   => $analysis,
        -epigenome  => $epigenome,
        -dbID       => 1000,
        -adaptor    => $issa,
        -experiment => $exp, );
}
qr/The feature type was undefined/, "Test constructor's feature_type exception";

throws_ok {
    my $set = Bio::EnsEMBL::Funcgen::Set->new(
        -name         => 'set_name',
        -feature_type => $feature_type,
        # -analysis     => $analysis,
        -epigenome  => $epigenome,
        -dbID       => 1000,
        -adaptor    => $issa,
        -experiment => $exp, );
}
qr/Set Analysis/, "Test constructor's analysis exception";

throws_ok {
    my $set =
        Bio::EnsEMBL::Funcgen::Set->new( -name         => 'set_name',
                                         -feature_type => $feature_type,
                                         -analysis     => $analysis,
                                         -epigenome    => 'not valid',
                                         -dbID         => 1000,
                                         -adaptor      => $issa,
                                         -experiment   => $exp, );
}
qr/Set Epigenome/, "Test constructor's epigenome exception";

throws_ok {
    my $set =
        Bio::EnsEMBL::Funcgen::Set->new( -name         => 'set_name',
                                         -feature_type => $feature_type,
                                         -analysis     => $analysis,
                                         -epigenome    => $epigenome,
                                         -dbID         => 1000,
                                         -adaptor      => $issa,
                                         -experiment   => 'not_valid', );
}
qr/Set Experiment/, "Test constructor's experiment exception";

# ------------
# Test getters
# ------------
is( $set->name,         'set_name',    "Test Set::name getter" );
is( $set->epigenome,    $epigenome,    "Test Set::epigenome getter" );
is( $set->feature_type, $feature_type, "Test Set::feature_type getter" );
is( $set->analysis,     $analysis,     "Test Set::analysis getter" );
is( $set->experiment_id, $set->{experiment_id},
    "Test Set::experiment_id getter" );

my $iss = $issa->fetch_by_name('SRR037563');
is( $iss->set_type, 'inputsub', "Test Set::set_type getter" );

is( $set->experiment, $exp, "Test Set::experiment getter" );

my $new_set =
    Bio::EnsEMBL::Funcgen::Set->new( -name          => 'set_name',
                                     -feature_type  => $feature_type,
                                     -analysis      => $analysis,
                                     -epigenome     => $epigenome,
                                     -dbID          => 1000,
                                     -adaptor       => $issa,
                                     -experiment_id => 4, );

is( $set->experiment, $exp,
    "Test Set::experiment getter (using experiment_id)" );

is( $set->source_label, 'ENCODE', "Test Set::source_label() getter" );

done_testing();
