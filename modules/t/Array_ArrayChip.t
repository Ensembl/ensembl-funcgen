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
use Bio::EnsEMBL::Funcgen::ArrayChip;
use Bio::EnsEMBL::Funcgen::Array;

use Test::More;
use Test::Exception;  # throws_ok
use Bio::EnsEMBL::Test::TestUtils qw( test_getter_setter debug );
use Bio::EnsEMBL::Test::MultiTestDB;

local $| = 1;

# switch on the debug prints
our $verbose = 0;

ok(1, 'Start up test');

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $db = $multi->get_DBAdaptor( "funcgen" );

isa_ok($db, 'Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor', 'Test database instatiated');


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

isa_ok($array, 'Bio::EnsEMBL::Funcgen::Array', 'Array constructor return type');

#
# What arrays should be able to do
# Temporarily set some new attributes

# These are really storable tests
#ok( test_getter_setter( $array, "dbID", 2 ), 'test_getter_setter dbID');
#ok( test_getter_setter( $array, "adaptor", undef ), 'test_getter_setter adaptor');

ok( test_getter_setter( $array, 'name', 'Test2' ), 'test_getter_setter Array::name');
ok( test_getter_setter( $array, 'format', 'TILED2' ), 'test_getter_setter Array::format');
ok( test_getter_setter( $array, 'vendor', 'ILLUMINA2' ), 'test_getter_setter Array::vendor');
ok( test_getter_setter( $array, 'class', 'ILLUMINA_WG2' ), 'test_getter_setter Array::classs');
ok( test_getter_setter( $array, 'description', 'TESTING2'), 'test_getter_setter Array::description');



#
#Create an ArrayChip
#

my $array_chip = Bio::EnsEMBL::Funcgen::ArrayChip->new 
  (-design_id   => 'Test',
   -array_id    => 34,  # HG-U133A
   -name        => 'Test',
   -adaptor     => $db->get_ArrayChipAdaptor); # Required for get_Array with -array_id

isa_ok($array_chip, 'Bio::EnsEMBL::Funcgen::ArrayChip', 'ArrayChip constructor return type (-array_id)');
# This currently throws as now adaptor has been set yet. This is only ever used in the ArrayChipAdaptor::_objs_from_sth for lazy loading

isa_ok($array_chip->get_Array, 'Bio::EnsEMBL::Funcgen::Array', 'ArrayChip::get_Array return type (-array_id)');

#Now use different params
$array_chip = Bio::EnsEMBL::Funcgen::ArrayChip->new 
  (-design_id   => 'Test_design_id',
   -array       => $array,
   -name        => 'Test_name');


#
# What array_chips should be able to do
#

#ok( test_getter_setter( $array_chip, "dbID", 2 ));
#ok( test_getter_setter( $array_chip, "adaptor", undef ));
isa_ok($array_chip, 'Bio::EnsEMBL::Funcgen::ArrayChip', 'ArrayChip constructor return type (-array)');
isa_ok( $array_chip->get_Array, 'Bio::EnsEMBL::Funcgen::Array', 'ArrayChip::get_Array return type (-array)');
is($array_chip->name, 'Test_name', 'ArrayChip::name');
is($array_chip->design_id, 'Test_design_id', 'ArrayChip::design_id');


#Now test some of the Array methods which require an ArrayChip
#Need to have adaptor here as the Array always queries the DB to
#Make sure we don't get any duplication of array chips


$array = $db->get_ArrayAdaptor->fetch_by_name_vendor('HG-U133A','AFFY');
#Only have 1 ArrayChip for this array
($array_chip) = @{$array->get_ArrayChips};
isa_ok($array_chip, 'Bio::EnsEMBL::Funcgen::ArrayChip', 'Array::get_ArrayChips return type');
is( scalar(@{$array->get_array_chip_ids}), 1, 'Array::get_array_chip_ids expected number');
my ($design_id) =  @{$array->get_design_ids};
is($design_id, 'HG-U133A', 'Array::get_design_ids expected ID');

$array_chip = $array->get_ArrayChip_by_design_id($design_id);
isa_ok($array_chip, 'Bio::EnsEMBL::Funcgen::ArrayChip', 'Array::get_ArrayChip_by_design_id return type');


$array_chip->{design_id} = 'TEST';  # design_id is only a getter
$array->add_ArrayChip($array_chip);
$array_chip = $array->get_ArrayChip_by_design_id('TEST');
isa_ok($array_chip, 'Bio::EnsEMBL::Funcgen::ArrayChip', 'Array::add_ArrayChip & get_ArrayChip_by_design_id');
is(scalar(@{$array->get_array_chip_ids}), 2, 'Array::add_ArrayChip & get_array_chip_ids returns expect amount');

done_testing();




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
