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
ok(scalar($sth->rows) == 17);


# now hide the reg_feat table i.e. make an empty version of it
$ens_test->hide("funcgen","regulatory_feature");
$sth->execute;
ok($sth->rows == 0);


# restore the reg feat table
$ens_test->restore();
$sth->execute;
ok(scalar($sth->rows) == 17);


# now save the reg feat table i.e. make a copy of it
$ens_test->save("funcgen","regulatory_feature");
$sth->execute;
ok(scalar($sth->rows) == 17);


# delete 5 reg feats from the db
$sth = $dba->dbc->prepare("delete from regulatory_feature where regulatory_feature_id <= 1617267 and regulatory_feature_id >=1617258");
$sth->execute;

$sth = $dba->dbc->prepare("select * from regulatory_feature");
$sth->execute;


ok(scalar($sth->rows) == 12);


# check to see whether the restore works again
$ens_test->restore();
$sth->execute;
ok(scalar($sth->rows) == 17);


$sth->finish;


