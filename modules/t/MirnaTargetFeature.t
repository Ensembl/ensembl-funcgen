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
BEGIN { use_ok('Bio::EnsEMBL::Funcgen::MirnaTargetFeature'); }

# ------------------------------
# Setup test database connection
# ------------------------------
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db    = $multi->get_DBAdaptor("funcgen");
my $fsa   = $db->get_adaptor("featureset");
my $fta   = $db->get_adaptor("featuretype");
my $aa    = $db->get_adaptor("analysis");
my $sa    = $db->get_adaptor("slice");
my $mtfa  = $db->get_adaptor("mirnatargetfeature");

# ----------------
# Test constructor
# ----------------
my $fs       = $fsa->fetch_by_name('A549_CTCF_ENCODE_Broad_SWEmbl_R0005_IDR');
my $ft       = $fta->fetch_by_name('H3K4me3');
my $analysis = $aa->fetch_by_logic_name('ChIP-Seq');
my $slice    = $sa->fetch_by_region( 'chromosome', 'X' );
my $start    = 1000;
my $end      = 2000;
my $strand   = 1;
my $display_label = 'test_label';

my $mtf = Bio::EnsEMBL::Funcgen::MirnaTargetFeature->new(
    -FEATURE_SET            => $fs,
    -FEATURE_TYPE           => $ft,
    -SLICE                  => $slice,
    -START                  => $start,
    -END                    => $end,
    -STRAND                 => $strand,
    -DISPLAY_LABEL          => $display_label,
    -ACCESSION              => 'test_accession',
    -ANALYSIS               => $analysis,
    -INTERDB_STABLE_ID      => 100,
    -EVIDENCE               => 'Computational',
    -METHOD                 => 'Microarray',
    -SUPPORTING_INFORMATION => 'info',

);

isa_ok( $mtf,
        'Bio::EnsEMBL::Funcgen::MirnaTargetFeature',
        'MirnaTargetFeature constructor return type' );

throws_ok {
    my $mtf = Bio::EnsEMBL::Funcgen::MirnaTargetFeature->new(
        -FEATURE_SET   => $fs,
        -FEATURE_TYPE  => $ft,
        -SLICE         => $slice,
        -START         => $start,
        -END           => $end,
        -STRAND        => $strand,
        -DISPLAY_LABEL => $display_label,
        # -ACCESSION     => 'test_accession',
        -ANALYSIS          => $analysis,
        -INTERDB_STABLE_ID => 100,

    );
}
qr/Mandatory parameter -accession not defined/,
    "Test constructor's exception for mandatory accession parameter";

# ------------
# Test getters
# ------------
is( $mtf->interdb_stable_id(),
    100, 'Test MirnaTargetFeature::interdb_stable_id() getter' );
is( $mtf->display_label(), $display_label,
    'Test MirnaTargetFeature::display_label() getter' );

my $new_mtf = Bio::EnsEMBL::Funcgen::MirnaTargetFeature->new(
    -FEATURE_SET  => $fs,
    -FEATURE_TYPE => $ft,
    -SLICE        => $slice,
    -START        => $start,
    -END          => $end,
    -STRAND       => $strand,
    # -DISPLAY_LABEL => $display_label,
    -ACCESSION         => 'test_accession',
    -ANALYSIS          => $analysis,
    -INTERDB_STABLE_ID => 100,
    -ADAPTOR           => $mtfa,

);

is( $new_mtf->display_label(),
    'CTCF - A549H3K4me3',
    'Test MirnaTargetFeature::display_label() getter without predefined label'
);

is( $mtf->accession(), 'test_accession',
    'Test MirnaTargetFeature::accession() getter' );

is( $mtf->evidence(), 'Computational',
    'Test MirnaTargetFeature::evidence() getter' );

is( $mtf->method(), 'Microarray',
    'Test MirnaTargetFeature::method() getter' );

is( $mtf->supporting_information(),
    'info', 'Test MirnaTargetFeature::supporting_information() getter' );

done_testing();
