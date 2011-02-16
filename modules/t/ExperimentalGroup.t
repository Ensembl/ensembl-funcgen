use strict;
use warnings;

use Bio::EnsEMBL::Funcgen::ExperimentalGroup;
use Bio::EnsEMBL::Funcgen::Experiment;

use Bio::EnsEMBL::Test::TestUtils qw( test_getter_setter debug );
use Bio::EnsEMBL::Test::MultiTestDB;

BEGIN { $| = 1;
	use Test;
	plan tests => 18;
}

# switch on the debug prints
our $verbose = 0;

debug( "Startup test" );
#1
ok(1);


my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();

my $db   = $multi->get_DBAdaptor( "funcgen" );
my $ega  = $db->get_ExperimentalGroupAdaptor;
my $ea =  $db->get_ExperimentAdaptor;

debug( "Test database instantiated" );
#2
ok( $db );

# The 'efg; should be there from the beginning
my $efg_group = $ega->fetch_by_name('efg');
#3
ok ( $efg_group );

#4
ok( $efg_group->name eq 'efg' );

my $group = Bio::EnsEMBL::Funcgen::ExperimentalGroup->new(
							   -name         => 'ebi_test',
							   -location     => 'location',
							   -contact      => 'contact',
							   -description  => 'Just a test group',
							 );
   

#5
ok( ref( $group ) && $group->isa( "Bio::EnsEMBL::Funcgen::ExperimentalGroup" ));

#6-11 Test
ok( test_getter_setter( $group, "dbID", 2 ));
ok( test_getter_setter( $group, "adaptor", undef ));
ok( test_getter_setter( $group, "name", "test name" ));
ok( test_getter_setter( $group, "location", "test location" ));
ok( test_getter_setter( $group, "contact", "test contact" ));
ok( test_getter_setter( $group, "description", "test description"));

#Prepare table for store tests....
$multi->hide('funcgen', 'experimental_group', 'experiment');
$ega->store($group);

$group = $ega->fetch_by_name("ebi_test");
#12-15 Test
ok( $group->name eq 'ebi_test' );
ok( $group->location eq 'location' );
ok( $group->contact eq 'contact' );
ok( $group->description =~ /^Just/ );

my $exp = Bio::EnsEMBL::Funcgen::Experiment->new
      (
       -NAME => 'test_experiment',
       -EXPERIMENTAL_GROUP => $group,
      );

$ea->store($exp);

$exp = $ea->fetch_by_name('test_experiment');

#16-18 Test
ok( $exp->name eq 'test_experiment' );
ok( $exp->experimental_group->name eq 'ebi_test' );
ok( $exp->experimental_group->description =~ /^Just/ );

#Restore table
$multi->restore();

1;
