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
use Bio::EnsEMBL::Funcgen::DBSQL::ProbeAdaptor;

# ------------------------------
# Setup test database connection
# ------------------------------
my $multi      = Bio::EnsEMBL::Test::MultiTestDB->new();
my $func_db    = $multi->get_DBAdaptor("funcgen");
my $probe_adaptor = Bio::EnsEMBL::Funcgen::DBSQL::ProbeAdaptor->new($func_db);

my $probe_set_array_name          = 'HG-U133A';
my $probe_name_on_probe_set_array = '552:605;';
my $probe_set_name                = '214727_at';

my $regular_array_name = 'WholeGenome_4x44k_v1';
my $regular_probe_name = 'A_24_P917810';

my $probe_adaptor = Bio::EnsEMBL::Registry->get_adaptor('homo_sapiens', 'Funcgen', 'Probe');

isa_ok( $probe_adaptor, 'Bio::EnsEMBL::Funcgen::DBSQL::ProbeAdaptor', 'Got a ProbeAdaptor');

my $probe_from_probe_set_fetched_the_old_way 
  = $probe_adaptor->fetch_by_array_probe_probeset_name(
  
    $probe_set_array_name, 
    $probe_name_on_probe_set_array, 
    $probe_set_name
  );

isa_ok( $probe_from_probe_set_fetched_the_old_way, 'Bio::EnsEMBL::Funcgen::Probe', 'Got a Probe');

is( 
  $probe_from_probe_set_fetched_the_old_way->name, 
  $probe_name_on_probe_set_array, 
  'name' 
);
is( 
  $probe_from_probe_set_fetched_the_old_way->array_chip->get_Array->name, 
  $probe_set_array_name, 
  'array name'
);
is( 
  $probe_from_probe_set_fetched_the_old_way->probe_set->name, 
  $probe_set_name, 
  'probe set'
);

my $probe_fetched_the_old_way 
  = $probe_adaptor->fetch_by_array_probe_probeset_name(
    $regular_array_name, 
    $regular_probe_name
  );

isa_ok( $probe_fetched_the_old_way, 'Bio::EnsEMBL::Funcgen::Probe', 'Got a Probe');

is( 
  $probe_fetched_the_old_way->name, 
  $regular_probe_name, 
  'name'
);

is( 
  $probe_fetched_the_old_way->array_chip->get_Array->name, 
  $regular_array_name, 
  'array name'
);

my $probe_fetched_the_new_way = $probe_adaptor->fetch_by_array_probe_name(
  $regular_array_name, 
  $regular_probe_name
);

isa_ok( 
  $probe_fetched_the_new_way, 
  'Bio::EnsEMBL::Funcgen::Probe', 
  'Got a Probe'
);

ok( 
  $probe_fetched_the_new_way->dbID == $probe_fetched_the_old_way->dbID, 
  'samesame'
);

my $probe_from_probe_set_fetched_the_new_way 
  = $probe_adaptor->fetch_by_array_probe_set_probe_name(

      $probe_set_array_name, 
      $probe_set_name,
      $probe_name_on_probe_set_array, 
  );

ok( 
  $probe_from_probe_set_fetched_the_new_way->dbID == $probe_from_probe_set_fetched_the_old_way->dbID, 
  'samesame'
);

done_testing();
