use strict;
use warnings;

use Bio::EnsEMBL::Funcgen::ExperimentalGroup;
use Bio::EnsEMBL::Funcgen::Experiment;

use Bio::EnsEMBL::Test::TestUtils qw( test_getter_setter debug );
use Bio::EnsEMBL::Test::MultiTestDB;

BEGIN { $| = 1;
	use Test;
	plan tests => 37;
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
							   -url          => 'http://www.ebi.ac.uk/',
							   -description  => 'Just a test group',
							   -is_project	 => 0
							 );
   

#5
ok( ref( $group ) && $group->isa( "Bio::EnsEMBL::Funcgen::ExperimentalGroup" ));

#6-13 Test
ok( test_getter_setter( $group, "dbID", 2 ));
ok( test_getter_setter( $group, "adaptor", undef ));
ok( test_getter_setter( $group, "name", "test name" ));
ok( test_getter_setter( $group, "location", "test location" ));
ok( test_getter_setter( $group, "contact", "test contact" ));
ok( test_getter_setter( $group, "url", "test url" ));
ok( test_getter_setter( $group, "description", "test description"));
ok( test_getter_setter( $group, "is_project", 1));

#Prepare table for store tests....
$multi->hide('funcgen', 'experimental_group', 'experiment');
$ega->store($group);

$group = $ega->fetch_by_name("ebi_test");
#14-19 Test
ok( $group->name eq 'ebi_test' );
ok( $group->location eq 'location' );
ok( $group->contact eq 'contact' );
ok( $group->url =~ /www.ebi/ );
ok( $group->description =~ /^Just/ );
ok( ! $group->is_project );

my $exp = Bio::EnsEMBL::Funcgen::Experiment->new
      (
       -NAME                => 'test_experiment',
       -EXPERIMENTAL_GROUP  => $group,
       -DATE                => "2011-01-01",
       -PRIMARY_DESIGN_TYPE => "test design",
       -ARCHIVE_ID	    => "GSEXXX",
       -DATA_URL	    => "http://",
       -DESCRIPTION         => "test description",     
      );


#20
ok( ref( $exp ) && $exp->isa( "Bio::EnsEMBL::Funcgen::Experiment" ));

#21-28 Test
ok( test_getter_setter( $exp, "dbID", 2 ));
ok( test_getter_setter( $exp, "adaptor", undef ));
ok( test_getter_setter( $exp, "experimental_group", $efg_group ));
ok( test_getter_setter( $exp, "date", "1999-01-01" ));
ok( test_getter_setter( $exp, "primary_design_type", "another test" ));
ok( test_getter_setter( $exp, "archive_id", "test id" ));
ok( test_getter_setter( $exp, "data_url", "test url" ));
ok( test_getter_setter( $exp, "description", "another description"));

$ea->store($exp);

$exp = $ea->fetch_by_name('test_experiment');

#29-37 Test
ok( $exp->name eq 'test_experiment' );
ok( $exp->experimental_group->name eq 'ebi_test' );
ok( $exp->experimental_group->description =~ /^Just/ );
ok( $exp->experimental_group->url =~ /www.ebi/ );
ok( $exp->date eq '2011-01-01' );
ok( $exp->primary_design_type eq 'test design' );
ok( $exp->archive_id eq 'GSEXXX' );
ok( $exp->data_url eq 'http://' );
ok( $exp->description eq 'test description' );

#Restore table
$multi->restore();

1;
