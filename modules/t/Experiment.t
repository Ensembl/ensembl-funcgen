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
use diagnostics;

use Data::Dumper qw( Dumper );
use Test::More;
use Test::Exception;    # throws_ok
use Bio::EnsEMBL::Test::TestUtils qw( test_getter_setter debug );

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::ExperimentalGroup;

# ---------------
# Module compiles
# ---------------
BEGIN { use_ok('Bio::EnsEMBL::Funcgen::Experiment'); }

# ------------------------------
# Setup test database connection
# ------------------------------
my $multi       = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db          = $multi->get_DBAdaptor("funcgen");
my $exp_adaptor = $db->get_ExperimentAdaptor;

# ----------------
# Test constructor
# ----------------
my $group = Bio::EnsEMBL::Funcgen::ExperimentalGroup->new(
    -name        => 'ebi_test',
    -location    => 'location',
    -contact     => 'contact',
    -url         => 'http://www.ebi.ac.uk/',
    -description => 'Just a test group',
    -is_project  => 0
);

my $exp = Bio::EnsEMBL::Funcgen::Experiment->new(
    -NAME                => 'test_experiment',
    -EXPERIMENTAL_GROUP  => $group,
    -DATE                => '2015-09-01',
    -PRIMARY_DESIGN_TYPE => 'test design',
    -ARCHIVE_ID          => 'GSEXXX',
    -DISPLAY_URL         => 'http://',
    -IS_CONTROL          => 1,
    -DESCRIPTION         => 'test description',
    -EPIGENOME           => $db->get_EpigenomeAdaptor->fetch_by_name('CD4'),
    -FEATURE_TYPE => $db->get_FeatureTypeAdaptor->fetch_by_name('CTCF'),
);

isa_ok(
    $exp,
    'Bio::EnsEMBL::Funcgen::Experiment',
    'Experiment constructor return type'
);

throws_ok {
    my $exp = Bio::EnsEMBL::Funcgen::Experiment->new(

        # -NAME                => 'test_experiment',
        -EXPERIMENTAL_GROUP  => $group,
        -DATE                => '2015-09-01',
        -PRIMARY_DESIGN_TYPE => 'test design',
        -ARCHIVE_ID          => 'GSEXXX',
        -DISPLAY_URL         => 'http://',
        -DESCRIPTION         => 'test description',
        -CELL_TYPE    => $db->get_CellTypeAdaptor->fetch_by_name('CD4'),
        -FEATURE_TYPE => $db->get_FeatureTypeAdaptor->fetch_by_name('CTCF'),
    );
}
qr/You must provide a name parameter/,
    'Test that a name parameter is provided to the constructor';

# ------------
# Test getters
# ------------
is_deeply(
    $exp->epigenome,
    $db->get_EpigenomeAdaptor->fetch_by_name('CD4'),
    'Test epigenome() getter'
);

is_deeply(
    $exp->feature_type,
    $db->get_FeatureTypeAdaptor->fetch_by_name('CTCF'),
    'Test feature_type() getter'
);

is( $exp->name, 'test_experiment', 'Test name() getter' );

is_deeply( $exp->experimental_group(),
    $group, 'Test experimental_group() getter' );

is_deeply( $exp->get_ExperimentalGroup(),
    $group, 'Test get_ExperimentalGroup() getter' );

is( $exp->description, 'test description', 'Test description() getter' );

# is( $exp->primary_design_type, 'test design',
#     'Test primary_design_type() getter' );

is( $exp->archive_id, 'GSEXXX', 'Test archive_id() getter' );

# is( $exp->display_url, 'http://', 'Test display_url() getter' );

# ----------------------
# Test getters - setters
# ----------------------
# ok( test_getter_setter( $exp, 'mage_xml', 'test xml' ),
#     'test_getter_setter Experiment::mage_xml'
# );

# ok( test_getter_setter( $exp, 'mage_xml_id', 5 ),
#     'test_getter_setter Experiment::mage_xml_id'
# );

# ---------------------
# Test get_InputSubsets
# ---------------------
my $new_exp = $exp_adaptor->fetch_by_name('NHLF_H3K9ac_ENCODE_Broad');

my $expected_is = {};    # expected input subsets
foreach my $input_subset (
    @{ $db->get_InputSubsetAdaptor->fetch_all_by_Experiments( [$new_exp] ) } )
{
    $expected_is->{ $input_subset->dbID } = $input_subset;
}

my @sorted_got = sort @{ $new_exp->get_InputSubsets };
my @sorted_expected = sort @{ [ values %{$expected_is} ] };

is( @sorted_got, @sorted_expected, 'Test get_InputSubsets()' );

# ------------------------------------------------------------------
# Test get_ExperimentalChips, add_ExperimentalChip,
# get_ExperimentalChip_by_unique_id, get_ExperimentalChip_unique_ids
# ------------------------------------------------------------------

# ********************************************#
# The above subroutines can not be tested,    #
# experimental_chip db table is always empty. #
# ********************************************#

# ------------------
# Test source_info()
# ------------------
my $expected_source_info = [[
    'ENCODE',
    'http://genome.ucsc.edu/cgi-bin/hgFileSearch?db=hg19&tsName=Broad&tsDescr=&tsGroup=regulation&fsFileType=fastq&hgt_mdbVar1=cell&hgt_mdbVal1=NHLF&hgt_mdbVar2=antibody&hgt_mdbVal2=H3K9ac&hgfs_Search=search'
]];

# is_deeply($new_exp->source_info(), $expected_source_info, 'Test source_info()');

# ----------------------------------
# Test reset_relational_attributes()
# ----------------------------------
my $new_cell_type = $db->get_EpigenomeAdaptor->fetch_by_name('U2OS');

my $new_exp_group = Bio::EnsEMBL::Funcgen::ExperimentalGroup->new(
    -name        => 'new_ebi_test',
    -location    => 'location',
    -contact     => 'contact',
    -url         => 'http://www.ebi.ac.uk/',
    -description => 'Just another test group',
    -is_project  => 0
);

my $new_feature_type = $db->get_FeatureTypeAdaptor->fetch_by_name('DNase1');

my $params = {
    -EPIGENOME          => $new_cell_type,
    -EXPERIMENTAL_GROUP => $new_exp_group,
    -FEATURE_TYPE       => $new_feature_type
};

my $existing_adaptor = $new_exp->adaptor();
my $existing_dbID    = $new_exp->dbID();

$new_exp->reset_relational_attributes( $params, 1 );

is_deeply(
    [   $new_exp->cell_type,    $new_exp->experimental_group,
        $new_exp->feature_type, $new_exp->adaptor(),
        $new_exp->dbID,
    ],
    [   $new_cell_type,    $new_exp_group, $new_feature_type,
        $existing_adaptor, $existing_dbID
    ],
    'Test reset_relational_attributes()'
);

$new_exp->reset_relational_attributes($params);

is_deeply(
    [ $new_exp->adaptor(), $new_exp->dbID, ],
    [ undef,               undef ],
    'Test reset_relational_attributes(), adaptor and dbID are undefined as expected'
);

done_testing();
