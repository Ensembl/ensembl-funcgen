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

use Test;
use strict;

BEGIN { $| = 1; plan tests => 9 }

use Bio::EnsEMBL::Test::MultiTestDB;

ok(1);

# Database will be dropped when this
# object goes out of scope
my $ens_test = Bio::EnsEMBL::Test::MultiTestDB->new;

ok($ens_test);

my $dba = $ens_test->get_DBAdaptor("funcgen");

ok($dba);


my $sth = $dba->dbc->prepare("select * from regulatory_feature");
$sth->execute;
ok(scalar($sth->rows) == 1);


# now hide the reg_feat table i.e. make an empty version of it
$ens_test->hide("funcgen","regulatory_feature");
$sth->execute;
ok($sth->rows == 0);


# restore the reg feat table
$ens_test->restore();
$sth->execute;
ok(scalar($sth->rows) == 1);


# now save the reg feat table i.e. make a copy of it
$ens_test->save("funcgen","regulatory_feature");
$sth->execute;
ok(scalar($sth->rows) == 1);


# delete 5 reg feats from the db
$sth = $dba->dbc->prepare("delete from regulatory_feature where regulatory_feature_id = 687807");
$sth->execute;

$sth = $dba->dbc->prepare("select * from regulatory_feature");
$sth->execute;


ok(scalar($sth->rows) == 0);


# check to see whether the restore works again
$ens_test->restore();
$sth->execute;
ok(scalar($sth->rows) == 1);


$sth->finish;


