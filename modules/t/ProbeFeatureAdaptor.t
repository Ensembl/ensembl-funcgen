# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2019] EMBL-European Bioinformatics Institute
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
use Bio::EnsEMBL::Registry;
use Test::More;
use Test::Exception;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Test::MultiTestDB;

my $multi      = Bio::EnsEMBL::Test::MultiTestDB->new();
my $func_db    = $multi->get_DBAdaptor("funcgen");

my $probe_set_array_name          = 'HuGene-2_0-st-v1';
my $probe_set_array_chip_name     = 'HuGene-2_0-st-v1';

my $probe_name_on_probe_set_array = '2398055-16657436';
my $probe_set_name                = '16657436';

my $regular_array_name      = 'OneArray';
my $regular_array_chip_name = 'OneArray';
my $regular_probe_name      = 'PH_hs_0000014';

my $probe_adaptor         = Bio::EnsEMBL::Registry->get_adaptor( 'homo_sapiens', 'Funcgen', 'Probe'        );
my $probe_set_adaptor     = Bio::EnsEMBL::Registry->get_adaptor( 'homo_sapiens', 'Funcgen', 'ProbeSet'     );
my $probe_feature_adaptor = Bio::EnsEMBL::Registry->get_adaptor( 'homo_sapiens', 'Funcgen', 'ProbeFeature' );

isa_ok( $probe_adaptor,         'Bio::EnsEMBL::Funcgen::DBSQL::ProbeAdaptor',        'Got a ProbeAdaptor'        );
isa_ok( $probe_set_adaptor,     'Bio::EnsEMBL::Funcgen::DBSQL::ProbeSetAdaptor',     'Got a ProbeSetAdaptor'     );
isa_ok( $probe_feature_adaptor, 'Bio::EnsEMBL::Funcgen::DBSQL::ProbeFeatureAdaptor', 'Got a ProbeFeatureAdaptor' );

my $probe_set = $probe_set_adaptor->fetch_by_array_probe_set_name(
  $probe_set_array_name,
  $probe_set_name,
);

isa_ok( $probe_set, 'Bio::EnsEMBL::Funcgen::ProbeSet', 'Got a ProbeSet');

my $probe_features_1 
  = $probe_feature_adaptor
    ->fetch_all_by_ProbeSet($probe_set)
;
my $probe_features_2 
  = $probe_feature_adaptor
    ->fetch_all_by_probeset_name(
        $probe_set_name
      )
;
my $probe_features_3 
  = $probe_feature_adaptor
    ->fetch_all_by_array_chip_name_probeset_name(
        $probe_set_array_chip_name, 
        $probe_set_name
      )
;
my $probe_features_4
  = $probe_feature_adaptor
    ->fetch_all_by_array_name_probeset_name(
        $probe_set_array_name, 
        $probe_set_name
      )
;

is(scalar @$probe_features_1, 57, 'Got expected number of features');
is(scalar @$probe_features_2, 57, 'Got expected number of features');
is(scalar @$probe_features_3, 57, 'Got expected number of features');
is(scalar @$probe_features_4, 57, 'Got expected number of features');

my $probe_features_5
  = $probe_feature_adaptor
    ->fetch_all_by_array_name_probe_name(
        $regular_array_name, 
        $regular_probe_name
      )
;
my $probe_features_6
  = $probe_feature_adaptor
    ->fetch_all_by_array_chip_name_probe_name(
        $regular_array_chip_name, 
        $regular_probe_name
      )
;

is(scalar @$probe_features_5, 3, 'Got expected number of features');
is(scalar @$probe_features_6, 3, 'Got expected number of features');

done_testing();
