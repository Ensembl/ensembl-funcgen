# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

use Test::More;
use Test::Exception;
use Bio::EnsEMBL::Funcgen::ExperimentalGroup;
use Bio::EnsEMBL::Funcgen::Experiment;
#use Bio::EnsEMBL::Funcgen::Utils::EFGUtils   qw( get_date );
use Bio::EnsEMBL::Test::TestUtils            qw( test_getter_setter debug );
use Bio::EnsEMBL::Test::MultiTestDB;

local $| = 1;
our $verbose = 0;  # switch on the debug prints

ok(1, 'Startup test');
my $multi     = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db        = $multi->get_DBAdaptor( "funcgen" );
my $ega       = $db->get_ExperimentalGroupAdaptor;
my $ea        =  $db->get_ExperimentAdaptor;
my $efg_group = $ega->fetch_by_name('Publication');

isa_ok($efg_group, 'Bio::EnsEMBL::Funcgen::ExperimentalGroup', 'ExperimentalGroupAdaptor::fetch_by_name return type');
is($efg_group->name, 'Publication', 'ExperimentalGroup::name');

my $group = Bio::EnsEMBL::Funcgen::ExperimentalGroup->new(
							   -name         => 'ebi_test',
							   -location     => 'location',
							   -contact      => 'contact',
							   -url          => 'http://www.ebi.ac.uk/',
							   -description  => 'Just a test group',
							   -is_project	 => 0
							 );
   
isa_ok($group, 'Bio::EnsEMBL::Funcgen::ExperimentalGroup', 'ExperimentalGroup::new return type');

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

#my $date = get_date('date');
my $exp = Bio::EnsEMBL::Funcgen::Experiment->new
 (-NAME                => 'test_experiment',
  -EXPERIMENTAL_GROUP  => $group,
  #-DATE                => $date,
  -PRIMARY_DESIGN_TYPE => 'test design',
  -ARCHIVE_ID	         => 'GSEXXX',
  -DISPLAY_URL         => 'http://',
  -DESCRIPTION         => 'test description',
  -CELL_TYPE           => $db->get_CellTypeAdaptor->fetch_by_name('CD4'),
  -FEATURE_TYPE        => $db->get_FeatureTypeAdaptor->fetch_by_name('CTCF'),     
 );



isa_ok($exp, 'Bio::EnsEMBL::Funcgen::Experiment', 'Experiment::new return type');
# These are only getters
is($exp->experimental_group, $group, 'Experiment::experimental_group');
#is($exp->date, $date, 'Experiment::date');
is($exp->primary_design_type, 'test design', 'Experiment::primary_design_type');
is($exp->archive_id, 'GSEXXX', 'Experiment::archive_id');
is($exp->display_url, 'http://', 'Experiment::display_url');
is($exp->description, 'test description', 'Experiment::description');

$ea->store($exp);

# This currently fails as ExperimentalGroup is not stored
$exp = $ea->fetch_by_name('test_experiment');

#29-37 Test
ok( $exp->name eq 'test_experiment' );
ok( $exp->experimental_group->name eq 'ebi_test' );
ok( $exp->experimental_group->description =~ /^Just/ );
ok( $exp->experimental_group->url =~ /www.ebi/ );
#is($exp->date, $date, 'Experiment::date after store/fetch');
ok( $exp->primary_design_type eq 'test design' );
ok( $exp->archive_id eq 'GSEXXX' );
ok( $exp->display_url eq 'http://' );
ok( $exp->description eq 'test description' );

#Restore table
$multi->restore();

done_testing();

1;
