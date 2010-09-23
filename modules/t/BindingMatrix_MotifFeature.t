use strict;
use warnings;

use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Funcgen::BindingMatrix;
use Bio::EnsEMBL::Funcgen::MotifFeature;
use Bio::EnsEMBL::Variation::VariationFeature;

use Bio::EnsEMBL::Test::TestUtils qw( test_getter_setter debug );
use Bio::EnsEMBL::Test::MultiTestDB;

BEGIN { $| = 1;
	use Test;
	plan tests => 46;
}

# switch on the debug prints
our $verbose = 0;

debug( "Startup test" );
#1
ok(1);

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $db   = $multi->get_DBAdaptor( "funcgen" );
my $bma  = $db->get_BindingMatrixAdaptor;
my $mf_a = $db->get_MotifFeatureAdaptor;

debug( "Test database instantiated" );
#2
ok( $db );


#
# Create BindingMatrix
#

my $ctcf_mat = ">CTCF_Test  
A [ 87	167	281	56	8	744	40	107	851	5	333	54	12	56	104	372	82	117	402  ]
C [ 291	145	49	800	903	13	528	433	11	0	3	12	0	8	733	13	482	322	181]
G [ 76	414	449	21	0	65	334	48	32	903	566	504	890	775	5	507	307	73	266 ] 
T [459	187	134	36	2	91	11	324	18	3	9	341	8	71	67	17	37	396	59 ] ";

my $ftype_a = $db->get_FeatureTypeAdaptor;
my $ftype   = $ftype_a->fetch_by_name('CTCF');

#3
ok ($ftype);

my $matrix = Bio::EnsEMBL::Funcgen::BindingMatrix->new(
													   -name         => 'CTCF',
													   -type         => 'Jaspar',
													   -feature_type => $ftype,
													   -description  => 'Jaspar Matrix',
													   -frequencies  => $ctcf_mat
													  );
   
#frequencies should be turned into _matrix 'blob'
#Maintain both methods:
#   frequencies for object creation?
#   matrix for utilisation. 


#4
ok( ref( $matrix ) && $matrix->isa( "Bio::EnsEMBL::Funcgen::BindingMatrix" ));

#
# What arrays should be able to do
# Temporarily set some new attributes

#5-9 Test
ok( test_getter_setter( $matrix, "dbID", 2 ));
ok( test_getter_setter( $matrix, "adaptor", undef ));
ok( test_getter_setter( $matrix, "name", "CTCF_TEST" ));
ok( test_getter_setter( $matrix, "type", "Inferred")); 
ok( test_getter_setter( $matrix, "description", "test description"));
# End of test 10


#Now hide tables for storage test
$multi->hide('funcgen', 'binding_matrix', 'associated_feature_type', 'motif_feature', 'associated_motif_feature');
$bma->store($matrix);
($matrix) = @{$bma->fetch_all_by_name("CTCF")}; #We know there is only one stored
#10
ok ( ref($matrix) && $matrix->isa('Bio::EnsEMBL::Funcgen::BindingMatrix') );

#Now test stored values
#11-14
#Test 12
ok( $matrix->name eq 'CTCF' );
ok( $matrix->type eq 'Jaspar' );
ok( $matrix->feature_type->name eq 'CTCF' );
ok( $matrix->description eq 'Jaspar Matrix' );
#End of test 15

#Need to write tests for these
#warn Data::Dumper::Dumper $matrix->_weights;
#warn $matrix->frequencies;

#numeric tests do not work with decimals?
#15-20
ok ( $matrix->_max_bind eq '18.100244783938' );
ok ( $matrix->relative_affinity("TGGCCACCAGGGGGCGCTA") == 1);
ok ( $matrix->relative_affinity("TGGCCACGAGGGGGCGCTA") eq '0.972088876164933');
ok ( $matrix->relative_affinity("TGGCCACCAGGGGGCGCCA") eq '0.997373533384315');
ok ( $matrix->relative_affinity("TGGCCACCAGGGGGCACTA") eq '0.99606869981799');
ok ( $matrix->relative_affinity("TGGCCACCAGGGAGCGCTA") eq '0.94541278085573');
#End of test 21

$matrix = Bio::EnsEMBL::Funcgen::BindingMatrix->new(
													-name  => 'CTCF',
													-type => 'Inferred',
													-description => 'Nkx3-2 Jaspar Matrix',
													-feature_type => $ftype,

												   );



$matrix->frequencies("A  [ 4  1 13 24  0  0  6  4  9 ]
C  [ 7  4  1  0  0  0  0  6  7 ]
G  [ 4  5  7  0 24  0 18 12  5 ]
T  [ 9 14  3  0  0 24  0  2  3 ]");

#Test  22 & 23
#21-23
ok ( $matrix->relative_affinity("TGGCCACCA") eq '0.209170634168983' );
ok ( $matrix->relative_affinity("CTCAGTGGA",1) eq '0.0655146380336576' );
ok( $matrix->length == 9);

#Now store second BindingMatrix
$bma->store($matrix);

#Test fetch methods
my @bms = @{$bma->fetch_all_by_name('CTCF')};

#24-25
ok (scalar (@bms) == 2);
ok (scalar (@{$bma->list_dbIDs}) == 2);

@bms = @{$bma->fetch_all_by_name('CTCF', 'inferred')};
#26
ok (scalar (@bms) == 1);

@bms = @{$bma->fetch_all_by_FeatureType($ftype)};
#27
ok (scalar (@bms) == 2);

#Now test stored values
#28-31
ok( $matrix->name eq 'CTCF' );
ok( $matrix->type eq 'Inferred' );
ok( $matrix->feature_type->name eq 'CTCF' );
ok( $matrix->description eq 'Nkx3-2 Jaspar Matrix' );


#Grab some random AFs for association
my $slice = $db->get_SliceAdaptor->fetch_by_region('chromosome', 1, 0, 2000000);
my ($af1, $af2)   = @{$db->get_AnnotatedFeatureAdaptor->fetch_all_by_Slice_FeatureType($slice, $ftype)};#CTCF

#Should really store these to avoid this
die('Test slice does not return enough AnnotatedFeatures, please change test set up') if (! ($af1 && $af2));


#Now let's test the MotifFeatures

my $mf = Bio::EnsEMBL::Funcgen::MotifFeature->new(
												  -binding_matrix => $matrix,
												  -score          => 10.1,
												  -start          => $af1->start,
												  -end            => $af1->end,
												  -slice          => $af1->slice,
												  -strand         => 1,
												 );
#32-35
ok( ref($mf) && $mf->isa('Bio::EnsEMBL::Funcgen::MotifFeature') );
#Now test the getter/setters
ok( test_getter_setter( $mf, 'display_label', 'test') );
ok( test_getter_setter( $mf, 'adaptor',       undef ) );
ok( test_getter_setter( $mf, 'score',            1.5) );

   
#Now simple store and fetch tests
($mf) = @{$mf_a->store($mf)};

#36
ok( ref($mf) && $mf->isa('Bio::EnsEMBL::Funcgen::MotifFeature') );

$mf = $mf_a->fetch_by_dbID($mf->dbID);
#37-40
ok( ref($mf) && $mf->isa('Bio::EnsEMBL::Funcgen::MotifFeature') );
# Test the existing and default attrs
ok( $mf->binding_matrix->dbID ==  $matrix->dbID );
ok( $mf->score                ==           10.1 );
ok( $mf->display_label  =~ /CTCF/ );


#Now store another MF and make associations

my $mf2 = Bio::EnsEMBL::Funcgen::MotifFeature->new(
												   -binding_matrix => $matrix,
												   -score          => 10.2,
												   -start          => $af2->start+1000,
												   -end            => $af2->start+1000+$matrix->length-1,
												   -slice          => $af2->slice,
												   -strand         => -1,
												 );

($mf2) = @{$mf_a->store($mf2)};


$mf  = $mf_a->store_associated_AnnotatedFeature($mf,  $af1);
$mf2 = $mf_a->store_associated_AnnotatedFeature($mf2, $af2);

#41
ok( $mf->associated_annotated_features->[0]->dbID == $af1->dbID );


#Test other fetch methods

my @mfs = @{$mf_a->fetch_all_by_AnnotatedFeature($af1)};
#42
ok( (scalar(@mfs) == 1) && ($mfs[0]->associated_annotated_features->[0]->dbID == $af1->dbID) );

ok($mf2->seq eq 'TGGACGCGG');
ok(length($mf2->seq) ==  $matrix->length);
ok($mf2->binding_matrix->relative_affinity($mf2->seq) eq '0.393642051097555');
my $mf_slice = $db->get_SliceAdaptor()->fetch_by_region('toplevel',$mf2->seq_region_name,$mf2->start,$mf2->end,$mf2->strand);
ok($mf2->seq eq $mf_slice->seq);

#Add a variation object
my $new_vf = Bio::EnsEMBL::Variation::VariationFeature->new(
  -start => 4,
  -end => 4,
  -slice => $mf->slice,       # the variation must be attached to a slice
  -allele_string => 'T/G',    # the first allele should be the reference allele
  -strand => -1,
  #-map_weight => 1,
  #-adaptor => $vfa,           # we must attach a variation feature adaptor
  -variation_name => 'newSNP',
);

ok($mf2->infer_variation_consequence($new_vf) eq '-18.0977682274379');

my %fsets;

foreach my $af($af1, $af2){
  $fsets{$af->feature_set->dbID} = $af->feature_set;
}

@mfs = @{$mf_a->fetch_all_by_Slice_FeatureSets($slice, [values %fsets])};
#43
ok( scalar(@mfs) == 2 );


@mfs = @{$mf_a->fetch_all_by_Slice_BindingMatrix($slice, $matrix)};
#44
ok( scalar(@mfs) == 2 );


my $ctype_id = $af1->cell_type->dbID;
my $ctype_count = 0;

foreach my $af($af1, $af2){
  $ctype_count ++ if $af->cell_type->dbID == $ctype_id;
}

@mfs = @{$mf_a->fetch_all_by_Slice_CellType($slice, $af1->cell_type)};
#45
ok( scalar(@mfs) == $ctype_count );


#list dbIDs
#46
ok( scalar(@{$mf_a->list_dbIDs}) == 2 );

#Restore table
$multi->restore();



1;
