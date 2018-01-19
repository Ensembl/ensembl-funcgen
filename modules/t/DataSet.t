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

use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis;
use Test::More;
use Test::Exception;               # throws_ok
use Bio::EnsEMBL::Test::TestUtils  qw( test_getter_setter debug );
use Data::Dumper                   qw( Dumper );
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Test::MultiTestDB;

# switch on the debug prints
our $verbose = 0;

# ---------------
# Module compiles
# ---------------
BEGIN { use_ok('Bio::EnsEMBL::Funcgen::DataSet'); }

# ------------------------------
# Setup test database connection
# ------------------------------
my $multi                = Bio::EnsEMBL::Test::MultiTestDB->new();
my $efgdba               = $multi->get_DBAdaptor("funcgen");
my $data_set_adaptor     = $efgdba->get_adaptor("dataset");
my $feature_set_adaptor  = $efgdba->get_adaptor("featureset");
my $result_set_adaptor   = $efgdba->get_adaptor("resultset");
my $analysis_adaptor     = $efgdba->get_adaptor("analysis");
my $epigenome_adaptor    = $efgdba->get_adaptor("epigenome");
my $feature_type_adaptor = $efgdba->get_adaptor("featuretype");

# ----------------
# Test constructor
# ----------------
my $result_set = $result_set_adaptor->fetch_by_name(
    'HeLa-S3_CTCF_ENCODE_Broad_bwa_samse');
my $feature_set = $feature_set_adaptor->fetch_by_name(
    'HeLa-S3_CTCF_ENCODE_Broad_SWEmbl_R0005_IDR');

my $supporting_sets=[$result_set];

my $data_set = Bio::EnsEMBL::Funcgen::DataSet->new(
    -NAME            => "test_name",
    -SUPPORTING_SETS => $supporting_sets,
    -FEATURESET      => $feature_set,
    -DISPLAYBLE      => 5,
);


isa_ok(
    $data_set,
    'Bio::EnsEMBL::Funcgen::DataSet',
    'DataSet constructor return type'
);

throws_ok {
    my $data_set = Bio::EnsEMBL::Funcgen::DataSet->new(

        # -NAME            => "test_name",
        -SUPPORTING_SETS => [$result_set],
        -FEATURESET      => $feature_set,
        -DISPLAYBLE      => 5,
    );
}
qr/Must defined a DataSet -name/, 'Test that a name is provided';

# ----------------------------------
# Test product_FeatureSet subroutine
# ----------------------------------
my $new_feature_set = $feature_set_adaptor->fetch_by_name(
    'HeLa-S3_CTCF_ENCODE_Uta_SWEmbl_R0005_IDR');
ok( test_getter_setter( $data_set, 'product_FeatureSet', $new_feature_set ),
    'test_getter_setter DataSet::product_FeatureSet' );

my $not_a_feature_set = $result_set;

throws_ok {
    $data_set->product_FeatureSet($not_a_feature_set);
}
qr/Need to pass a valid Bio::EnsEMBL::Funcgen::FeatureSet/,
    'Test that a FeatureSet object is provided to product_FeatureSet()';

throws_ok {
    $data_set->product_FeatureSet($feature_set);
}
qr/The main feature_set has already been set for this DataSet, maybe you want add_SupportingSets?/,
    'Test that the new FeatureSet object provided to product_FeatureSet() is not the same as the existing one';

# -----------------------------------------------
# Test get_supporting_sets_by_Analysis subroutine
# -----------------------------------------------
my $analysis = $analysis_adaptor->fetch_by_logic_name('bwa_samse');

is_deeply( $data_set->get_supporting_sets_by_Analysis($analysis),
    $supporting_sets, 'Test get_supporting_sets_by_Analysis() subroutine' );

my $not_an_analysis = $result_set;

throws_ok {
    $data_set->get_supporting_sets_by_Analysis($not_an_analysis);
}
qr/Need to pass a valid stored Bio::EnsEMBL::Analysis/,
    'Test that a Bio::EnsEMBL::Analysis object is provided to get_supporting_sets_by_Analysis()';

# -----------------------------------
# Test get_supporting_sets subroutine
# -----------------------------------
is_deeply( $data_set->get_supporting_sets(),
    $supporting_sets, 'Test get_supporting_sets() subroutine' );

throws_ok {
    $data_set->get_supporting_sets('invalid_set_type');
}
qr/You have specified an invalid supporting set type/,
    'Test that a valid $set_type is provided to get_supporting_sets()';

# -----------------------------------------------
# Test get_displayable_supporting_sets subroutine
# -----------------------------------------------
#is_deeply( $data_set->get_displayable_supporting_sets,
#    $supporting_sets, 'Test get_displayable_supporting_sets() subroutine' );

# --------------------------------------------------
# Test get_displayable_product_FeatureSet subroutine
# --------------------------------------------------
#is_deeply( $data_set->get_displayable_product_FeatureSet,
#    $new_feature_set,
#    'Test get_displayable_product_FeatureSet() subroutine' );

# -------------------------------------------------------
# Test name(), epigenome() and feature_type() subroutines
# -------------------------------------------------------
is( $data_set->name(), 'test_name', 'Test name() subroutine' );

my $expected_epigenome = $epigenome_adaptor->fetch_by_name('HeLa-S3');
is_deeply( $data_set->epigenome(), $expected_epigenome,
    'Test cell_type() subroutine' );

my $expected_feature_type = $feature_type_adaptor->fetch_by_name('CTCF');
is_deeply( $data_set->feature_type(),
    $expected_feature_type, 'Test feature_type() subroutine' );

# -------------------------------
# Test display_label() subroutine
# -------------------------------
is( $data_set->display_label(),
    'CTCF - HeLa-S3 Enriched Sites',
    'Test display_label() subroutine'
);


# # This test uses the following comment conventions
# # START method_name testing
# # COMPLETED method_name testing

# #Just grab a few data sets to work with
# #DISPLAYABLE ensure we should have an input_set defined
# #$data_set_adaptor->fetch_all;
# my $dset   = $data_set_adaptor->fetch_by_name('RegulatoryFeatures:MultiCell');
# my $dset_2 = $data_set_adaptor->fetch_by_name('RegulatoryFeatures:NHEK');

# if(! (defined $dset && defined $dset_2)){
#   throw('Failed to fetch 2 DataSets to test, please update DataSet.t');  
# }


# # START testing compare_to

# my $diffs = '';
# my %diffs = %{$dset->compare_to($dset)};
# $diffs = "\n".Dumper(\%diffs) if %diffs;
# ok(! %diffs, 'DataSet::compare_to self default no diffs'.$diffs);

# # redefine methods
# %diffs = %{$dset->compare_to($dset, undef,
#   [qw(name get_all_states)] ,
#   [qw(product_FeatureSet get_supporting_sets)]  
#   )};
# $diffs = "\n".Dumper(\%diffs) if %diffs;
# ok(! %diffs, 'DataSet::compare_to self redefined methods no diffs'.$diffs);

# eval{#Invalid redefined scalar method
#  $dset->compare_to($dset, undef, 
#                    [qw(product_FeatureSet get_supporting_sets)] ); 
# };
# ok($@, 'DataSet::compare_to scalar methods redefined (invalid)');


# eval{#Invalid redefined object method
#  $dset->compare_to($dset, undef, undef,
#                    [qw(name get_all_states)] ); 
# };
# ok($@, 'DataSet::compare_to object methods redefined (invalid)');

# # COMPLETED testing compare_to

# ## START testing reset_relational_attributes                                             

# #use eq directly on object refs instead of compare_to


# #clone ResultSet, so we can change some attrs
# my $clone_dset = bless({%{$dset}}, ref($dset));
# my $alt_fset   = $dset_2->product_FeatureSet;
# my $alt_ssets  = $dset_2->get_supporting_sets;

# my %relational_params = 
#   (
#    -feature_set     => $alt_fset,
#    -supporting_sets => $alt_ssets,
#   );


# $clone_dset->reset_relational_attributes(\%relational_params, 'no_db_reset');


# #eq comparisons of obj in scalar context compares mem refs
# ok( ($clone_dset->product_FeatureSet eq $alt_fset),
#    'DataSet::reset_relational_attributes reset FeatureSet');

# #can't compare arrayrefs directly, so use compare_to here
# %diffs = %{$clone_dset->compare_to($dset_2, undef, [],
#                                       ['get_supporting_sets'])};
# $diffs = "\n".Dumper(\%diffs) if %diffs;
# ok(! $diffs,
#    'DataSet::reset_relational_attributes reset supporting sets'.$diffs);


# ok((defined $clone_dset->dbID && 
#    defined $clone_dset->adaptor), 
#    'DataSet::reset_relational_attributes no_db_reset');
    
# $clone_dset->reset_relational_attributes(\%relational_params);
# ok(! (defined $clone_dset->dbID || 
#       defined $clone_dset->adaptor), 
#       'DataSet::reset_relational_attributes with dbID/adaptor reset');



# eval { $clone_dset->reset_relational_attributes(
#         {  
#          -feature_set     => $alt_fset,
#         });
# };
# ok($@, 'DataSet::reset_relational_attributes no -feature_set error');

# eval { $clone_dset->reset_relational_attributes(
#          {  
#           -supporting_sets => $alt_ssets,
#          })
# };
# ok($@, 'DataSet::reset_relational_attributes no -supporting_sets error');



#todo _set_Sets_and_types

done_testing();


1;

