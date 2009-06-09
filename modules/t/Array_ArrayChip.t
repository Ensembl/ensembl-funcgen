use strict;
use warnings;

use Bio::EnsEMBL::CoordSystem;
use Bio::EnsEMBL::Slice;
use Bio::EnsEMBL::Funcgen::ArrayChip;
use Bio::EnsEMBL::Funcgen::Array;

use Bio::EnsEMBL::Test::TestUtils qw( test_getter_setter debug );
use Bio::EnsEMBL::Test::MultiTestDB;

BEGIN { $| = 1;
	use Test;
	plan tests => 23;
}

# switch on the debug prints
our $verbose = 0;

debug( "Startup test" );
ok(1);

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $db = $multi->get_DBAdaptor( "funcgen" );

debug( "Test database instatiated" );
ok( $db );


#
# Create Array
#

my $array = Bio::EnsEMBL::Funcgen::Array->new 
  (
   -NAME        => 'Test',
   -FORMAT      => 'Tiled',
   -SIZE        => '1',
   -VENDOR      => 'VENDOR',
   -TYPE        => 'OLIGO',
   -DESCRIPTION => 'TESTING',
   -CLASS       => 'VENDOR_FORMAT',#e.g. AFFY_UTR, ILLUMINA_WG
  );

ok( ref( $array ) && $array->isa( "Bio::EnsEMBL::Funcgen::Array" ));

#
# What arrays should be able to do
# Temporarily set some new attributes

ok( test_getter_setter( $array, "dbID", 2 ));
ok( test_getter_setter( $array, "adaptor", undef ));
ok( test_getter_setter( $array, "name", "Test2" ));
ok( test_getter_setter( $array, "format", 'TILED2' ));
ok( test_getter_setter( $array, "vendor", 'ILLUMINA2' ));
ok( test_getter_setter( $array, "class", 'ILLUMINA_WG2' ));
ok( test_getter_setter( $array, "description", 'TESTING2'));
#10 tests


#
#Create an ArrayChip
#

#What about testing for failures?
#eval and test $@?

my $array_chip = Bio::EnsEMBL::Funcgen::ArrayChip->new 
  (
   -design_id   => 'Test',
   -array_id    => 1,
   -name        => 'Test',

  );

ok( ref( $array_chip ) && $array_chip->isa( "Bio::EnsEMBL::Funcgen::ArrayChip" ));

#Now use different params
$array_chip = Bio::EnsEMBL::Funcgen::ArrayChip->new 
  (
   -design_id   => 'Test',
   -array       => $array,
   -name        => 'Test',
  );

ok( ref( $array_chip ) && $array_chip->isa( "Bio::EnsEMBL::Funcgen::ArrayChip" ));

#
# What array_chips should be able to do
#

ok( test_getter_setter( $array_chip, "dbID", 2 ));
ok( test_getter_setter( $array_chip, "adaptor", undef ));
ok( test_getter_setter( $array_chip, "name", "Test2" ));
ok( test_getter_setter( $array_chip, "design_id", "Test2" ));
ok( ref($array_chip->get_Array) eq 'Bio::EnsEMBL::Funcgen::Array');
#17 tests


#Now test some of the Array methods which require an ArrayChip
#Need to have adaptor here as the Array always queries the DB to
#Make sure we don't get any duplication of array chips


$array = $db->get_ArrayAdaptor->fetch_by_name_vendor('HG-U133A','AFFY');
#Only have 1 ArrayChip for this array
($array_chip) = @{$array->get_ArrayChips};
ok( ref($array_chip) eq 'Bio::EnsEMBL::Funcgen::ArrayChip');
ok( scalar(@{$array->get_array_chip_ids}) == 1);
my ($design_id) =  @{$array->get_design_ids};
ok($design_id eq 'HG-U133A');

$array_chip = $array->get_ArrayChip_by_design_id($design_id);
ok( ref( $array_chip ) && $array_chip->isa( "Bio::EnsEMBL::Funcgen::ArrayChip" ));

#
$array_chip->design_id('TEST');
$array->add_ArrayChip($array_chip);
$array_chip = $array->get_ArrayChip_by_design_id('TEST');
ok( ref( $array_chip ) && $array_chip->isa( "Bio::EnsEMBL::Funcgen::ArrayChip" ));
ok( scalar(@{$array->get_array_chip_ids}) == 2);


#Num test 23





__END__

# possibly do the following with database connection
#my $probes = $affy_array->get_all_AffyProbes();

#
# Create a Probe
#
# Several ways to create

my $probe = Bio::EnsEMBL::Funcgen::Probe->new
  (
   -array => $array,
   -name => "123-145",
   -probeset => "affy_probeset"
  );

ok( ref( $probe ) && $probe->isa( "Bio::EnsEMBL::Funcgen::Probe" ));

my $merge_probe = Bio::EnsEMBL::AffyProbe->new
(
  -arraynames => [ "Affy-1", "Affy-2", "Affy-3" ],
  -probenames => [ "123-145", "23,24,56", "someplace" ],
  -probeset => "affy_probeset"
);

ok( ref( $merge_probe ) && $merge_probe->isa( "Bio::EnsEMBL::AffyProbe" ));




#
# What probes need to be able to do
#

ok( test_getter_setter( $merge_probe, "dbID", 1 ));
ok( test_getter_setter( $merge_probe, "adaptor", bless( {}, "Bio::EnsEMBL::DBSQL::BaseAdaptor" )));
ok( test_getter_setter( $merge_probe, "probeset", "Affy_probeset" ));


my $arrays = $affy_probe->get_all_AffyArrays();
ok( ref( $arrays ) eq "ARRAY" );
ok( $arrays->[0] && $arrays->[0]->isa( "Bio::EnsEMBL::AffyArray" ));


# expected to construct full probenames from Chipname, setname and probename
my $full_probenames = $merge_probe->get_all_complete_names();
ok( ref( $full_probenames ) eq 'ARRAY' );
ok( $full_probenames->[0] && ( $full_probenames->[0] =~ /affy_probeset/ ));
 

my $full_probename = $merge_probe->get_complete_name( "Affy-1" );
ok( $full_probename =~ /affy_probeset/ && $full_probename =~ /Affy-1/ );

my $probenames = $merge_probe->get_all_probenames();
my $probename = $merge_probe->get_probename( $affy_array->name() );


#
# When we implement storing this should be implemented as well
# for starters we are not creating these objects, but load the database with probes ...

# $affy_probe->add_array_name( $affy_array, $probename );


#
# Create an affy feature
#

my $coord_system = Bio::EnsEMBL::CoordSystem->new
   ( -name => "chromosome",
     -version => '',
     -rank => 1 );

my $slice = Bio::EnsEMBL::Slice->new
   ( 
     -seq_region_name => '1',
     -coord_system => $coord_system,
     -start => 1,
     -end => 50,
     -seq_region_length => 200_000_000
);
     
my $affy_feature = Bio::EnsEMBL::AffyFeature->new
    (
     -probe => $affy_probe,
     -mismatch_count => 1,
     -slice => $slice,
     -start => 1,
     -end => 10,
     -strand => -1,
     -probeset => "affy_probeset"	
);

ok( ref( $affy_feature ) && $affy_feature->isa( "Bio::EnsEMBL::AffyFeature"));

1;
