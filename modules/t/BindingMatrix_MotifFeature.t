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

use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Funcgen::BindingMatrix;
use Bio::EnsEMBL::Funcgen::MotifFeature;
use Bio::EnsEMBL::Variation::VariationFeature;

use Bio::EnsEMBL::Test::TestUtils qw( test_getter_setter debug );
use Bio::EnsEMBL::Test::MultiTestDB;
use Test::More;
use Test::Exception;  # throws_ok

# Test DB slice is 13:32888400-32974305
# switch on the debug prints
our $verbose = 0;

# TODO
# Add can_ok tests for all methods?

my $skip = 0;
ok(1, 'Start up');

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db    = $multi->get_DBAdaptor( 'funcgen' );
ok( $db, 'Test database instantiated');

my $bma  = $db->get_BindingMatrixAdaptor;
my $mf_a = $db->get_MotifFeatureAdaptor;


# Create BindingMatrix

my $freqs = ">CTCF_Test
A [ 87  167 281 56  8 744 40  107 851 5 333 54  12  56  104 372 82  117 402  ]
C [ 291 145 49  800 903 13  528 433 11  0 3 12  0 8 733 13  482 322 181]
G [ 76  414 449 21  0 65  334 48  32  903 566 504 890 775 5 507 307 73  266 ]
T [459  187 134 36  2 91  11  324 18  3 9 341 8 71  67  17  37  396 59 ] ";

my $ftype_a = $db->get_FeatureTypeAdaptor;
my $ftype_name = 'CTCF';
my $ftype   = $ftype_a->fetch_by_name($ftype_name);

my $bm_name = 'MA0139.1';
my $lname = 'JASPAR_v5_0';
my $bm_desc  = 'Jaspar Matrix';
my $analysis_a = $db->get_AnalysisAdaptor;
my $analysis   = $analysis_a->fetch_by_logic_name($lname);



my $matrix = Bio::EnsEMBL::Funcgen::BindingMatrix->new(
                             -name         => $bm_name,
                             -analysis     => $analysis,
                             -feature_type => $ftype,
                             -description  => $bm_desc,
                             -frequencies  => $freqs,
                             -threshold    => 0.1,
                            );
my $m_ctcf = $matrix;
   
# TODO Need to add test for 2d array passed to frequncies


isa_ok($matrix, 'Bio::EnsEMBL::Funcgen::BindingMatrix', 'BindingMatrix new');
ok($matrix->threshold == 0.1,                           'threshold stored from new');
ok($matrix->length == 19,                               'Matrix length from string');
ok(test_getter_setter($matrix, "threshold", 0.9),       'Set/get threshold attribute');
ok($matrix->is_position_informative(5),                 'Position 5 is informative');
ok(! $matrix->is_position_informative(19),              'Position 9 is not informative');

throws_ok(sub{ $matrix->is_position_informative(100) }, 
          qr/.*Position\([0-9]+\) is out of bounds.*/s, 
          'BindingMatrix is_position_informative out of bounds');

#Now hide tables for storage test
$multi->hide('funcgen', 'binding_matrix', 'associated_feature_type', 'motif_feature', 'associated_motif_feature');
$bma->store($matrix);
($matrix) = @{$bma->fetch_all_by_name($bm_name)}; #We know there is only one stored

$skip = ! isa_ok($matrix, 'Bio::EnsEMBL::Funcgen::BindingMatrix', 'Fetched stored BindingMatrix');
#ok (! $skip, 'Fetched stored BindingMatrix');

SKIP:{  # Test stored values
  $skip && skip "Failed to retrieve stored BindingMatrix", 6; 
  ok( $matrix->dbID =~ /^[0-9]+$/, 'Stored dbID returned');
  isa_ok($matrix->adaptor, 'Bio::EnsEMBL::Funcgen::DBSQL::BindingMatrixAdaptor', 
    'BindingMatrixAdaptor set');
  is($matrix->name, $bm_name, 'Stored name returned');
  is($matrix->description, $bm_desc, 'Stored description returned' );
  is($matrix->analysis->logic_name, $lname, 'Stored logic_name returned');
  is($matrix->feature_type->name, $ftype_name, 'Stored FeatureType name returned' );
}

$matrix = $m_ctcf; # in case store failed
#Need to write tests for these
#warn Data::Dumper::Dumper $matrix->_weights;
#warn $matrix->frequencies;

#numeric tests do not work with decimals?
#TODO change this to sprintf 3 decimal places! i.e. float not double
is( $matrix->max_affinity, '18.0919082398347', 'max_affinity');
is( $matrix->relative_affinity("TGGCCACCAGGGGGCGCTA"), 1, 'relative_affinity - consensus');
is( $matrix->relative_affinity("TGGCCACGAGGGGGCGCTA"), '0.972088876164933', 'relative_affinity - 8th C->G');
is( $matrix->relative_affinity("TGGCCACCAGGGGGCGCCA"), '0.997373533384315', 'relative_affinity - 18thT -> C');
is( $matrix->relative_affinity("TGGCCACCAGGGGGCACTA"), '0.99606869981799', 'relative_affinity - 16th G -> A');
is( $matrix->relative_affinity("TGGCCACCAGGGAGCGCTA"), '0.94541278085573', 'relative_affinity - 13th G -> A');

# New analysis so we can fetch 2 with the same name form diferent analyses
$analysis = $analysis_a->fetch_by_logic_name('Regulatory_Build');
$freqs = [
  [ qw(4  1 13 24  0  0  6  4  9) ],  # A
  [ qw(7  4  1  0  0  0  0  6  7) ],  # C
  [ qw(4  5  7  0 24  0 18 12  5) ],  # G
  [ qw(9 14  3  0  0 24  0  2  3) ]];  # T
$matrix = Bio::EnsEMBL::Funcgen::BindingMatrix->new(
                          -name         => $bm_name, #this is  a bit odd
                          -analysis     => $analysis,
                          -description  => 'Regulatory_Build Matrix',
                          -feature_type => $ftype,
                          -frequencies  => $freqs,
                           );

ok( ! defined $matrix->threshold, 'Threshold undef');
$matrix->threshold(0.81);

is($matrix->relative_affinity("TGGCCACCA"), '0.209170634168983', 'Relative affinity test (log)');
is($matrix->relative_affinity("CTCAGTGGA", 1), '0.0655146380336576', 'Relative affinity test (linear)');
is($matrix->relative_affinity("NTCAGTGGA", 1), undef, 'Relative affinity test seq contains N');
#Now store second BindingMatrix
$bma->store($matrix);

#Test fetch methods
my @bms = @{$bma->fetch_all_by_name($bm_name)};
is(scalar(@bms), 2, "Fetched duplicate $bm_name matrices by name");
#ok (scalar (@{$bma->_list_dbIDs}) == 2);  # Is this necessarily true? Also a BaseAdaptor test?

@bms = @{$bma->fetch_all_by_name($bm_name, $analysis)};
is(scalar(@bms), 1, "Fetch unique analysis specific $bm_name matrix");
@bms = @{$bma->fetch_all_by_FeatureType($ftype)};
is(scalar(@bms), 2, "Fetched duplicate $bm_name matrices by FeatureType");

is($matrix->length, 9, 'Corrected matrix length');

# Do we need these? Or is length test a good proxy
# as we have already tested these for above matrix
ok($matrix->is_position_informative(5), 'BindingMatrix position 5 is informative');
ok(! $matrix->is_position_informative(1), 'BindingMatrix position 1 is not informative');

# Grab some random AFs for association... avoid areas close to the edge!
# 6 1024000 1028940 This is the only region in the test DB and they are histones
# data is also on GRCh37, and DB probably needs patching?
# hmm, we need
# Does not help that the test DB has had the coord systems switched for some reason
# This is going to mess thing up if we copy things from the GRCh37 DB
# dbIDs are not goign to match and sequence dependant tests are going to be hard
# to handle i.e. won't be able to get a true peak for relative affinity calcs


my $slice = $db->get_SliceAdaptor->fetch_by_region('chromosome', 13, 32888400, 32975000);
my ($af1, $af2)   = @{$db->get_AnnotatedFeatureAdaptor->fetch_all_by_Slice_FeatureType($slice, $ftype)};#CTCF

# Test DB has 18 CTCF peaks


#SKIP:{
# ($af1 && $af2) || 
#   skip 'Test slice does not return enough AnnotatedFeatures, please change test set up', ;

#Now let's test the MotifFeatures

# this is creating a MotifFeature the same length as the af, not the bm!
# this should trigger errors is we call is_position_informative
# We don't expose relative_affinity via a mf wrapper

my $mf = Bio::EnsEMBL::Funcgen::MotifFeature->new(
                          -binding_matrix => $matrix,
                          -score          => 0.9,
                          -start          => $af1->start,
                          -end            => $af1->start + $matrix->length -1,
                          -slice          => $af1->slice,
                          -strand         => 1,
                          #no -display_label here as we test default below  
                         );
isa_ok($mf, 'Bio::EnsEMBL::Funcgen::MotifFeature', 'MotifFeature new');

#ok( test_getter_setter( $mf, 'adaptor',       undef ) ); # Not a MotifFeature test
($mf) = @{$mf_a->store($mf)};
$mf = $mf_a->fetch_by_dbID($mf->dbID);
isa_ok($mf, 'Bio::EnsEMBL::Funcgen::MotifFeature', 'MotifFeatureAdaptor::fetch_by_dbID');
# don't actually test dbID here

# Test the existing and default attrs
is($mf->display_label, 
   $mf->binding_matrix->feature_type->name.':'.$mf->binding_matrix->name, 
   'MotifFeature::display_label default');
# where is this set?


is($mf->binding_matrix->dbID,  $matrix->dbID, 'MotifFeature BindingMatrix dbID match');
is($mf->score, 0.9, 'MotifFeature::score set via new');

ok($mf->is_position_informative(5), 'MotifFeature position 5 is informative');

#warn "\n".$mf->binding_matrix->frequencies;

# This should be exactly the same as the BindingMatrix test above, why is this failing but that not?
# Value is the same?!!
ok(! $mf->is_position_informative(1), 'MotifFeature position 1 is not informative');


# Now store another downstream MF and make associations

my $mf2 = Bio::EnsEMBL::Funcgen::MotifFeature->new(
                           -binding_matrix => $matrix,
                           -score          => 0.2,
                           -start          => $af2->start + 1000,
                           -end            => $af2->end + 1000, #Intentio
                           -slice          => $af2->slice,
                           -strand         => -1, 
                           -display_label  => 'test',
                         );
($mf2) = @{$mf_a->store($mf2)};
is($mf2->display_label, 'test', 'MotifFeature::display_label set via new');
$mf  = $mf_a->store_associated_AnnotatedFeature($mf,  $af1);
$mf2 = $mf_a->store_associated_AnnotatedFeature($mf2, $af2);

is($mf->associated_annotated_features->[0]->dbID, $af1->dbID, 
  'Stored associated_annotated_feature dbID match');

#Test other fetch methods

my @mfs = @{$mf_a->fetch_all_by_AnnotatedFeature($af1)};

ok((scalar(@mfs) == 1) && ($mfs[0]->associated_annotated_features->[0]->dbID == $af1->dbID),
   'MotifFeatureAdaptor::fetch_all_by_AnnotatedFeature returns expected associated MotifFeature');


ok($mf2->is_position_informative(5), 'Position 5 revcomp is informative');
ok(! $mf2->is_position_informative(9), 'Position 9 revcomp is not informative');

#Add a variation object
my $vf = Bio::EnsEMBL::Variation::VariationFeature->new(
  -start         => 5,
  -end           => 5,
  -slice         => $mf->feature_Slice,       # the variation must be attached to a slice
  -allele_string => 'A/G',    # the first allele should be the reference allele, always wrt +ve strand
  -strand        => -1,  # what?
#  #-map_weight => 1,
#  #-adaptor => $vfa,           # we must attach a variation feature adaptor
  -variation_name => 'testSNP',
  #No attached variation here, so this may cause some odd failures
);


is($mf2->infer_variation_consequence($vf), 0, 'infer_variation_consequence out of bounds has no consequence');
# use stderr_is here?

# Using revcomp allels due to $vf strand
$vf->allele_string('T/C');
#cmp_ok($mf->infer_variation_consequence($vf), '>', 0, 'infer_variation_consequence 5 more weight');
is($mf->infer_variation_consequence($vf), '18.0977682274379', 'infer_variation_consequence 5 more weight');


# we also need to test with strand 1?
$vf = Bio::EnsEMBL::Variation::VariationFeature->new(
  -start         => 1,
  -end           => 1,
  -slice         => $mf->feature_Slice,     
  -allele_string => 'A/C',    
  -strand        => 1,  
  -variation_name => 'testSNP',
  #No attached variation here, so this may cause some odd failures
);

is($mf->infer_variation_consequence($vf), 
   '1.81184937270422', 
   'infer_variation_consequence 1 little consequence');

$vf->allele_string('A/G');
is($mf->infer_variation_consequence($vf), 0, 'infer_variation_consequence 1 equal weight');

# todo unsupported allele string

my %fsets;

foreach my $af($af1, $af2){
  $fsets{$af->feature_set->dbID} = $af->feature_set;
}

@mfs = @{$mf_a->fetch_all_by_Slice_FeatureSets($slice, [values %fsets])};
is(scalar(@mfs), 2, 'MotifFeatureAdaptor::fetch_all_by_Slice_FeatureSets expected number');

@mfs = @{$mf_a->fetch_all_by_Slice_BindingMatrix($slice, $matrix)};
is(scalar(@mfs), 2, 'MotifFeatureAdaptor::fetch_all_by_Slice_BindingMatrix expected number');

my $ctype_id = $af1->cell_type->dbID;
my $ctype_count = 0;

foreach my $af($af1, $af2){
  $ctype_count ++ if $af->cell_type->dbID == $ctype_id;
}

@mfs = @{$mf_a->fetch_all_by_Slice_CellType($slice, $af1->cell_type)};
is(scalar(@mfs), $ctype_count, 'MotifFeatureAdaptor::fetch_all_by_Slice_CellType expected number');

#Restore table
$multi->restore;

# N base test

TODO:{
  local $TODO = 'Seq contains N for MotifFeature::infer_variation_consequence tests';


  # Kinda hard to do this on MotifFeature as the test DB currently does not have any Ns!
  # mysql> select seq_region_id, length(sequence) l from dna where sequence like "%N%" order by l asc limit 10;
# +---------------+------+
# | seq_region_id | l    |
# +---------------+------+
# |    2006139307 |  100 |
# |        115151 |  358 |
# |    2006139387 |  395 |
# |    2006139385 |  792 |
# |    2006139318 |  940 |
# |    2006139379 | 1095 |
# |    2006139386 | 1140 |
# |         35438 | 1194 |
# |         35456 | 1408 |
# |         35454 | 1703 |
# +---------------+------+
# Need to insert one of these into the test DB, then create a MotifFeature over a region with an N
# Then do a relative_affinity call where the base change is on a normal base (should use neutral weight)
# Then do a relative_affinity where the base changed is the N base (should return undef)
# Does the VEP (MotifFeatureVariationAllele::motif_score_delta) even use infer_variation_consequence? No reference in
# ensembl-variation. So need to check code between motif_score_delta and infer_variation_consequence.

}

done_testing();

1;
