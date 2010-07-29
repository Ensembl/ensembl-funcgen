use strict;
use warnings;

use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Funcgen::BindingMatrix;
use Bio::EnsEMBL::Funcgen::MotifFeature;

use Bio::EnsEMBL::Test::TestUtils qw( test_getter_setter debug );
use Bio::EnsEMBL::Test::MultiTestDB;

BEGIN { $| = 1;
	use Test;
	plan tests => 34;
}

# switch on the debug prints
our $verbose = 0;

debug( "Startup test" );
ok(1);

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $db = $multi->get_DBAdaptor( "funcgen" );
my $bma = $db->get_BindingMatrixAdaptor(); 


debug( "Test database instantiated" );
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



ok( ref( $matrix ) && $matrix->isa( "Bio::EnsEMBL::Funcgen::BindingMatrix" ));

#
# What arrays should be able to do
# Temporarily set some new attributes

ok( test_getter_setter( $matrix, "dbID", 2 ));
ok( test_getter_setter( $matrix, "adaptor", undef ));
ok( test_getter_setter( $matrix, "name", "CTCF_TEST" ));
ok( test_getter_setter( $matrix, "type", "Inferred")); 
ok( test_getter_setter( $matrix, "description", "test description"));


#Now hide tables for storage test
$multi->hide('funcgen', 'binding_matrix', 'associated_feature_type', 'motif_feature');
$bma->store($matrix);
($matrix) = @{$bma->fetch_all_by_name("CTCF")}; #We know there is only one stored
ok ( ref($matrix) && $matrix->isa('Bio::EnsEMBL::Funcgen::BindingMatrix') );

#Now test stored values
ok( $matrix->name eq 'CTCF' );
ok( $matrix->type eq 'Jaspar' );
ok( $matrix->feature_type->name eq 'CTCF' );
ok( $matrix->description eq 'Jaspar Matrix' );


#Need to write tests for these
#warn Data::Dumper::Dumper $matrix->_weights;
#warn $matrix->frequencies;

#numeric tests do not work with decimals?
ok ( $matrix->_delta eq '72583253.4125366' );
ok ( $matrix->compare_to_optimal_site("TGGCCACCAGGGGGCGCTA") == 1);
ok ( $matrix->compare_to_optimal_site("TGGCCACGAGGGGGCGCTA") eq '0.11105980143154');
ok ( $matrix->compare_to_optimal_site("TGGCCACCAGGGGGCGCCA") eq '0.813178490280233');
ok ( $matrix->compare_to_optimal_site("TGGCCACCAGGGGGCACTA") eq '0.733780319463616');
ok ( $matrix->compare_to_optimal_site("TGGCCACCAGGGAGCGCTA") eq '0.0135939782046961');

#This is faked data to test duplicate names with different types.
my $assoc_ftype = $ftype_a->fetch_by_name('DNase1');


ok ($assoc_ftype);

$matrix = Bio::EnsEMBL::Funcgen::BindingMatrix->new(
													-name  => 'CTCF',
													-type => 'Inferred',
													-description => 'Nkx3-2 Jaspar Matrix',
													-feature_type => $ftype,
													-associated_feature_types => [$assoc_ftype],
                                                       );

$matrix->frequencies("A  [ 4  1 13 24  0  0  6  4  9 ]
C  [ 7  4  1  0  0  0  0  6  7 ]
G  [ 4  5  7  0 24  0 18 12  5 ]
T  [ 9 14  3  0  0 24  0  2  3 ]");

ok ( $matrix->compare_to_optimal_site("TGGCCACCA") eq '3.90079276657717e-11' );
ok ( $matrix->compare_to_optimal_site("CTCAGTGGA") eq '0.0655146380337265' );

#Now store second BindingMatrix
$bma->store($matrix);

#Test fetch methods
my @bms = @{$bma->fetch_all_by_name('CTCF')};

ok (scalar (@bms) == 2);
ok (scalar (@{$bma->list_dbIDs}) == 2);

@bms = @{$bma->fetch_all_by_name('CTCF', 'inferred')};
ok (scalar (@bms) == 1);

@bms = @{$bma->fetch_all_by_FeatureType($ftype)};
ok (scalar (@bms) == 2);

@bms = @{$bma->fetch_all_by_associated_FeatureType($assoc_ftype)};
ok (scalar (@bms) == 1);

$matrix = $bms[0];

ok( scalar(@{$matrix->associated_feature_types}) == 1);
ok( $matrix->associated_feature_types->[0]->name eq 'DNase1' );


#Now test stored values
ok( $matrix->name eq 'CTCF' );
ok( $matrix->type eq 'Inferred' );
ok( $matrix->feature_type->name eq 'CTCF' );
ok( $matrix->description eq 'Nkx3-2 Jaspar Matrix' );



#Restore table
$multi->restore();



1;
