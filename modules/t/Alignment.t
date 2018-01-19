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

use Test::More                     qw( no_plan );
use Data::Dumper                   qw( Dumper );
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Utils::Scalar    qw( check_ref );
use Bio::EnsEMBL::Funcgen::Alignment;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils; # warns_like test_getter_setter

use Test::Exception;  # throws_ok # This only work when Error objects are used

# ---------------
# Module compiles
# ---------------
BEGIN { use_ok('Bio::EnsEMBL::Funcgen::Alignment'); }

# ------------------------------
# Setup test database connection
# ------------------------------
my $multi  = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db     = $multi->get_DBAdaptor('funcgen');
isa_ok($db, 'Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor', 'Test database instantiated');


# This test uses the following comment conventions
# START method_name testing
# COMPLETED method_name testing


my $rsa     = $db->get_adaptor('resultset');
# my $slice_a = $db->dnadb->get_SliceAdaptor;
my $aa      = $db->get_adaptor("analysis");
my $exp_a   = $db->get_adaptor('experiment');
my $fta     = $db->get_adaptor("featuretype");
my $issa    = $db->get_adaptor("inputsubset");
my $epi_a   = $db->get_adaptor("epigenome");

# ----------------
# Test constructor
# ----------------
my $analysis       = $aa->fetch_by_logic_name('SWEmbl_R0005_IDR');
my $feature_type   = $fta->fetch_by_name('CTCF');
my $new_result_set = Bio::EnsEMBL::Funcgen::Alignment->new(
    -analysis      => $analysis,
    -feature_class => 'result',
    -feature_type  => $feature_type,
    -name          => 'new_result_set',
    -table_name    => 'input_subset',
    -replicate     => 3,
    -adaptor       => $rsa,
);

isa_ok(
    $new_result_set,
    'Bio::EnsEMBL::Funcgen::Alignment',
    'Alignment constructor result type'
);

throws_ok {
    my $new_result_set = Bio::EnsEMBL::Funcgen::Alignment->new(
        -analysis      => $analysis,
        -feature_class => 'result',
        -feature_type  => $feature_type,
        -name          => 'new_result_set',
        -table_name    => 'invalid table_name',
    );
}
qr/You need to pass a valid -table_name/,
    "Test exception throw for invalid table_name ";


my $iss = $issa->fetch_by_name('SRR037563');

throws_ok {
    my $new_result_set = Bio::EnsEMBL::Funcgen::Alignment->new(
        -analysis      => $analysis,
        -feature_class => 'result',
        -feature_type  => $feature_type,
        -name          => 'new_result_set',
        -table_id      => 50,
        -support       => [$iss],
    );
}
qr/Unsafe to specify -support and -table_id, please use -support/,
    "Test exception throw for support and table name parameter";



my $exp_name = 'H1ESC_Tcf12_ENCODE_Hudsonalpha';
my $exp      = $exp_a->fetch_by_name($exp_name);
 
# SKIP:{
#   if(! $exp){
#     skip "Could not fetch test Experiment:\t$exp_name\nThis must have been removed from the DB, please choose another Experiment or fix the test DB";
#   }
    
#   #eval { @rsets = @{$rsa->fetch_all_by_Experiment($exp)} };
#   #ok(! $@, "AlignmentAdaptor::fetch_all_by_Experiment\t$@");
#   # No point in evaling as we depend on @4rsets below
#   # Let's just let it fail here for now, until we can be bothered to skip appropriately?
#   #ok(@rsets = @{$rsa->fetch_all_by_Experiment($exp)}, 'AlignmentAdaptor::fetch_all_by_Experiment');
#   #This is mostly null and void now as it is testing ResultFeature stuff which has been striped out

#   my @rsets =  @{$rsa->fetch_all_by_Experiment($exp)};

#   is(scalar(@rsets), 1,
#    'AlignmentAdaptor::fetch_all_by_Experiment returns expected number of Alignments');

#   my $rset_exp       = $rsets[0]->experiment;
#   is($rset_exp->dbID, $exp->dbID, 'Alignment::get_Experiment returns expected Experiment');
    
 
#   #TODO
#   #Now let's get and check the file path for each window size

# }


# #Just grab a few result sets to work with
# my ($result_set, $result_set_2) = 
#   @{ $rsa->fetch_all_by_feature_class('dna_methylation', 
#                                       {status => 'DISPLAYABLE'}) };
# if(! (defined $result_set && defined $result_set_2)){
#   throw('Failed to fetch 2 Alignments to test, please update Alignment.t');  
# }

# # -----------------
# # Test compare_to()
# # -----------------
# # START testing compare_to
# my %diffs = %{$result_set->compare_to($result_set)};
# my $diffs = '';
# $diffs = "\n".Dumper(\%diffs) if %diffs;
# ok(! %diffs, 'Alignment::compare_to self defaults returns no diffs'.$diffs);

# # redefine methods
# %diffs = %{$result_set->compare_to($result_set, undef,
#   [qw(name table_name feature_class get_all_states)] ,
#   [qw(feature_type cell_type analysis get_support)]  
#   )};
# $diffs = "\n".Dumper(\%diffs) if %diffs;
# ok(! %diffs, 'Alignment::compare_to self redefined methods returns no diffs'.$diffs);


# # Unavailable redefined scalar method
# # This is actually a Storable::compare_scalar_methods test!
# throws_ok {$result_set->compare_to($result_set, undef, undef, [qw(unavailable_method)] )}
#   qr/cannot call method/s,
#   'Storable/Alignment::compare_to catches unavailable method';

# # Invalid redefined scalar method
# # This is actually a Storable::compare_scalar_methods test!
# throws_ok {$result_set->compare_to($result_set, undef, [qw(feature_type cell_type analysis get_support)] )}
#   qr/does not return a SCALAR value or an ARRAY or ARRAYREF of SCALAR values:/s,
#   'Storable/Alignment::compare_to catches invalid redefined scalar methods';

# # Invalid redefined object method
# # This is actually a Storable::_compare_method_return_types test!
# throws_ok {$result_set->compare_to($result_set, undef, undef, [qw(name table_name feature_class get_all_states)] )}
#   qr/ method does not return Storable\(s\)/s, 
#   'Storable/Alignment::compare_to catches invalid redefined storable methods';

# # COMPLETED testing compare_to

# # ----------------------------------
# # Test reset_relational_attributes()
# # ----------------------------------
# ## START testing reset_relational_attributes 

# # clone Alignment, so we can change some attrs
# my $clone_rset   = bless({%{$result_set}}, ref($result_set));
# my $alt_ftype    = $result_set->feature_type->adaptor->fetch_by_name('H3K4me3');
# my $alt_analysis = $result_set->analysis->adaptor->fetch_by_logic_name('AFFY_UTR_ProbeAlign');
# my $alt_ctype    = ($result_set_2->cell_type->name eq 'H1ESC') ?
#                     $result_set_2->cell_type->adaptor->fetch_by_name('GM06990') :
#                     $result_set_2->cell_type->adaptor->fetch_by_name('H1ESC');
# my $alt_exp      = $result_set_2->experiment;
# my $alt_support  = $result_set_2->get_support;
# my $orig_support = $result_set->get_support;#also lazyload now before we blow away the adaptor
# #todo check support is different!



# my %relational_params = 
#  (-support      => $alt_support,
#   -analysis     => $alt_analysis,
#   -feature_type => $alt_ftype,
#   -cell_type    => $alt_ctype,
#   -experiment   => $alt_exp     );

# $clone_rset->reset_relational_attributes(\%relational_params, 'no_db_reset');


# #eq comparisons of obj in scalar context compares mem refs
# ok( ($clone_rset->analysis eq $alt_analysis),
#    'Alignment::reset_relational_attributes reset analysis');


# #todo test support ID differ here

# ok( ($clone_rset->feature_type eq $alt_ftype),
#    'Alignment::reset_relational_attributes reset feature_type');

# ok( ($clone_rset->cell_type eq $alt_ctype),
#    'Alignment::reset_relational_attributes reset ctype');
   
# ok((defined $clone_rset->dbID && 
#    defined $clone_rset->adaptor), 
#    'Alignment::reset_relational_attributes no_db_reset');
    
# $clone_rset->reset_relational_attributes(\%relational_params);
# ok(! (defined $clone_rset->dbID || 
#       defined $clone_rset->adaptor), 
#       'Alignment::reset_relational_attributes with dbID/adaptor reset');


# eval { $clone_rset->reset_relational_attributes(
#          {  
#           -support      => $result_set->get_support,
#           -analysis     => $result_set->analysis,
#           -feature_type => $result_set->feature_type,
#          })
# };

# # We are not testing the error string here?

# ok($@, 'Alignment::reset_relational_attributes no -cell_type error');

# eval { $clone_rset->reset_relational_attributes(
#          {  
#           -support   => $result_set->get_support,
#           -analysis  => $result_set->analysis,
#           -cell_type => $result_set->cell_type,
#          })
# };
# ok($@, 'Alignment::reset_relational_attributes no -feature_type error');

# eval { $clone_rset->reset_relational_attributes(
#          {  
#           -support   => $result_set->get_support,
#           -cell_type => $result_set->cell_type,
#           -feature_type => $result_set->feature_type,
#          })
# };
# ok($@, 'Alignment::reset_relational_attributes no -analysis error');

# eval { $clone_rset->reset_relational_attributes(
#          {  
#           -cell_type => $result_set->cell_type,
#           -feature_type => $result_set->feature_type,
#           -analysis => $result_set->analysis,
#          })
# };
# ok($@, 'Alignment::reset_relational_attributes no -support error');



# # COMPLETED testing reset_relational_attributes

# # ------------------
# # Test add_support()
# # ------------------
# # START testing add_support
# #my @orig_support   = @{$orig_support};
# #my @orig_support_2 = @{$result_set_2->get_support};
# #push @orig_support, @orig_support_2;
# my @new_support    = @{$clone_rset->add_support($orig_support)};#This appends rather than resets!

# #let's assume if we have the same size, the we have been successful?
# ok( ( (scalar(@new_support)) == 1 &&
#       (scalar(@$orig_support))),
#     'Alignment::add_support returns expected size array');
    
# ok( ( ($new_support[0] eq $orig_support->[0])),
#     'Alignment::add_support returns expected object');

# #Duplicate support test has to be done on result_set_2 
# #to avoid borking the rest of the tests
# eval { $clone_rset->add_support($orig_support) };
# ok($@, 'Alignment::add_support caught duplicate support addition');  

# # ------------------
# # Test dbfile_path()
# # ------------------
# # START testing dbfile_path dbfile_data_root etc

# is($rsa->dbfile_data_root, '', 'AlignmentAdaptor default dbfile_data_root is null string');
# is($rsa->dbfile_data_root('/test/root'), '/test/root', 'AlignmentAdaptor::dbfile_data_root set');

# # This no longer updates dynamically
# is($result_set->dbfile_path, 
#   '/dna_methylation_feature/AG04449_5mC_ENCODE_Uw/wgEncodeHaibMethylRrbsAg04449UwstamgrowprotSitesfiltered_10.bb',
#   'Alignment::dbfile_path');

# # So we have to fetch it again now we have set the dbfile_data_root
# $result_set = $rsa->fetch_all_by_feature_class('dna_methylation', {status => 'DISPLAYABLE'})->[0];
# is($result_set->dbfile_path, 
#   '/test/root/dna_methylation_feature/AG04449_5mC_ENCODE_Uw/wgEncodeHaibMethylRrbsAg04449UwstamgrowprotSitesfiltered_10.bb',
#   'Alignment::dbfile_path');
# # hide table, store this result set and retest returned values

# ----------------------
# Test production_name()
# ----------------------
# TODO Without production names for epigenomes in the test db this is never going to pass
# is($new_result_set->production_name(),'', 'Test production_name()');

# -----------------------------
# Test _valid_feature_classes()
# -----------------------------
my @valid_feature_classes = $new_result_set->_valid_feature_classes();
my @expected_classes = ('result','dna_methylation','segmentation');
is(@valid_feature_classes, @expected_classes,'Test _valid_feature_classes()');

# ----------------------------------
# Test reset_relational_attributes()
# ----------------------------------
# my $epigenome        = $epi_a->fetch_by_name('H7ESC');
# my $reset_parameters = {
#     -analysis     => $analysis,
#     -feature_type => $feature_type,
#     -epigenome    => $epigenome
# };
# # my $reset_parameters = [$analysis,$feature_type, $epigenome];


# throws_ok {
#     $new_result_set->reset_relational_attributes($reset_parameters);
# }
# qr/You must pass a valid Bio::EnsEMBL::Funcgen::Experiment/,
#     "Test exception throw for reset_relational_attributes()";


# ----------------
# Test replicate()
# ----------------
# is($new_result_set->replicate(),3,'Test replicate()');

# ------------------
# Test add_support()
# ------------------
# TODO 

# --------------------
# Test display_label()
# --------------------
# TODO 

# -----------------
# Test table_name()
# -----------------
is($new_result_set->table_name(),'input_subset','Test table_name()');

# --------------------
# Test _add_table_id()
# --------------------
$new_result_set->_add_table_id(5, 1);
is($new_result_set->{table_ids}->{'5'},1, 'Test _add_table_id()' );

throws_ok {
    $new_result_set->_add_table_id( 5, 1 );
}
qr/You are attempting to redefine a result_set_input_id which is already defined/,
    "Test exception throw for already defined result_set_input_id";

# ----------------
# Test table_ids()
# ----------------
is_deeply($new_result_set->table_ids(), [5], 'Test table_ids()');

# ---------------------------
# Test result_set_input_ids()
# ---------------------------
is_deeply( $new_result_set->result_set_input_ids(),
    [1], 'Test result_set_input_ids()' );

# ---------------
# Test contains()
# ---------------
# TODO

# ------------------------------
# Test get_result_set_input_id()
# ------------------------------
is($new_result_set->get_result_set_input_id(5), 1, 'Test get_result_set_input_id()');
is($new_result_set->get_result_set_input_id(4), undef, 'Test get_result_set_input_id() again');

# ------------------
# Test get_support()
# ------------------
my $brand_new_result_set = Bio::EnsEMBL::Funcgen::Alignment->new(
        -analysis      => $analysis,
        -feature_class => 'result',
        -feature_type  => $feature_type,
        -name          => 'new_result_set',
        -support       => [$iss],
    );
is_deeply($brand_new_result_set->get_support(),[$iss], 'Test get_support()');
