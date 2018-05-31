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

use Test::More;
use Test::Exception;

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::RegulatoryBuild;
use Bio::EnsEMBL::Funcgen::FeatureType;
use Bio::EnsEMBL::Analysis;



my $multi   = Bio::EnsEMBL::Test::MultiTestDB->new();
my $func_db = $multi->get_DBAdaptor('funcgen');
my $core_db = $multi->get_DBAdaptor('core');

# ---------------
# Module compiles
# ---------------
BEGIN { use_ok('Bio::EnsEMBL::Funcgen::DBSQL::RegulatoryBuildAdaptor'); }

my $regulatoryBuild_adaptor         = $func_db->get_adaptor('RegulatoryBuild');
my $featureType_adaptor         = $func_db->get_adaptor('FeatureType');
my $analysis_adaptor         = $func_db->get_adaptor('Analysis');



# ----------------
# Test adaptor
# ----------------
isa_ok( $regulatoryBuild_adaptor, 'Bio::EnsEMBL::Funcgen::DBSQL::RegulatoryBuildAdaptor',
    'adaptor: RegulatoryBuildAdaptor' );

# ----------------
# Test fetch_by_name
# ----------------
my $regulatory_build= $regulatoryBuild_adaptor->fetch_by_name('Regulatory features');
isa_ok( $regulatory_build, 'Bio::EnsEMBL::Funcgen::RegulatoryBuild',
    'fetch_by_name: RegulatoryBuild' );

is($regulatory_build->name, 'Regulatory features', 'fetch_by_name');

$regulatory_build= $regulatoryBuild_adaptor->fetch_by_name('Non existant name');
is($regulatory_build, undef, 'fetch_by_name: Non existant name');

# ----------------
# Test fetch_current_regulatory_build
# ----------------
my $current_regulatory_build= $regulatoryBuild_adaptor->fetch_current_regulatory_build();
isa_ok( $current_regulatory_build, 'Bio::EnsEMBL::Funcgen::RegulatoryBuild',
    'fetch_current_regulatory_build: Current_RegulatoryBuild' );

$current_regulatory_build= $regulatoryBuild_adaptor->fetch_current_regulatory_build();
ok($current_regulatory_build->is_current(),'fetch_current_regulatory_build: is_current');

# ----------------
# Test store
# ----------------

my $regulatory_build_logic_name        = 'Regulatory_Build';
my $regulatory_build_feature_type_name = 'RegulatoryFeature';

my $regulatory_build_analysis = $analysis_adaptor->fetch_by_logic_name($regulatory_build_logic_name);
my $regulatory_build_feature_type = $featureType_adaptor->fetch_by_name($regulatory_build_feature_type_name);

my $new_regulatoryBuild = Bio::EnsEMBL::Funcgen::RegulatoryBuild->new(
    -name                         => 'new RegBuild Name',
    -feature_type              => $regulatory_build_feature_type,
    -analysis                  => $regulatory_build_analysis,
    -is_current                   => 1,
    );

$regulatoryBuild_adaptor->store($new_regulatoryBuild);

my $fetched_regulatoryBuild = $regulatoryBuild_adaptor->fetch_by_name('new RegBuild Name');
my $fetched_FeatureType = $fetched_regulatoryBuild->fetch_FeatureType;
my $fetched_Analysis =$fetched_regulatoryBuild->fetch_Analysis;

isa_ok( $fetched_regulatoryBuild, 'Bio::EnsEMBL::Funcgen::RegulatoryBuild', 'store: fetched RegulatoryBuild' );
is($fetched_regulatoryBuild->name, 'new RegBuild Name', "store: fetched RegulatoryBuild Name");
isa_ok( $fetched_FeatureType, 'Bio::EnsEMBL::Funcgen::FeatureType', 'store: fetched FeatureType' );
is( $fetched_FeatureType->name, $regulatory_build_feature_type_name, 'store: fetched FeatureType Name' );
isa_ok( $fetched_Analysis, 'Bio::EnsEMBL::Analysis', 'store: fetched Analysis' );
is( $fetched_Analysis->logic_name, $regulatory_build_logic_name, 'store: Analysis LogicName' );
is( $fetched_regulatoryBuild->is_current, 1, 'store: fetched RegulatoryBuild is_current' );

$regulatory_build = $regulatoryBuild_adaptor->fetch_by_name('Regulatory features');
is( $regulatory_build->is_current, 0, 'store: existing RegulatoryBuild is_current' );

# ----------------
# Test update
# ----------------

my $regulatory_build_dbID = $fetched_regulatoryBuild->dbID;
$fetched_regulatoryBuild->{dbID} = undef;
throws_ok { $regulatoryBuild_adaptor->update($fetched_regulatoryBuild) }
qr/The database id must be set/,
  'update: fetched RegulatoryBuild No dbID exception throw';


$fetched_regulatoryBuild->dbID($regulatory_build_dbID);
$fetched_regulatoryBuild->name('Updated Name');
$fetched_regulatoryBuild->version('Updated Version');
$fetched_regulatoryBuild->initial_release_date('18/05/2018');
$fetched_regulatoryBuild->last_annotation_update('19/05/2018');
$fetched_regulatoryBuild->feature_type_id(4);
$fetched_regulatoryBuild->analysis_id(5);
$fetched_regulatoryBuild->is_current(1);
$fetched_regulatoryBuild->sample_regulatory_feature_id(6);

$regulatoryBuild_adaptor->update($fetched_regulatoryBuild);

$fetched_regulatoryBuild = $regulatoryBuild_adaptor->fetch_by_name('Updated Name');

isa_ok( $fetched_regulatoryBuild, 'Bio::EnsEMBL::Funcgen::RegulatoryBuild', 'update: fetched RegulatoryBuild' );

is($fetched_regulatoryBuild->name, 'Updated Name', "update: fetched RegulatoryBuild Name");
is($fetched_regulatoryBuild->version, 'Updated Version', "update: fetched RegulatoryBuild Version");
is($fetched_regulatoryBuild->initial_release_date, '18/05/2018', "update: fetched RegulatoryBuild initial_release_date");
is($fetched_regulatoryBuild->last_annotation_update, '19/05/2018', "update: fetched RegulatoryBuild last_annotation_update");
is($fetched_regulatoryBuild->feature_type_id, 4, "update: fetched RegulatoryBuild feature_type_id");
is($fetched_regulatoryBuild->analysis_id, 5, "update: fetched RegulatoryBuild analysis_id");
is($fetched_regulatoryBuild->is_current, 1, "update: fetched RegulatoryBuild is_current");
is($fetched_regulatoryBuild->sample_regulatory_feature_id, 6, "update: fetched RegulatoryBuild sample_regulatory_feature_id");

$regulatory_build->is_current(1);
$regulatoryBuild_adaptor->update($regulatory_build);
$fetched_regulatoryBuild = $regulatoryBuild_adaptor->fetch_by_name('Updated Name');

is( $regulatory_build->is_current, 1, 'update: RegulatoryBuild is_current' );
is( $fetched_regulatoryBuild->is_current, 0, 'update: existing RegulatoryBuild is_current' );



done_testing();