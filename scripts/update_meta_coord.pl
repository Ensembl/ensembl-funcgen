#!/usr/local/ensembl/bin/perl -w

use strict;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Getopt::Long;

my $help = 0;
my ($pass, $dbpattern );
my $user = 'ensadmin';
my $port = 3306;
my $host = 'ens-genomics1';

my $usage = qq(
  $0 --host ens-staging --port 3306 --user ensadmin \\
    --pass XXX --dbpattern core

  [--help] displays this menu.

This script will dump the current meta_coord table in the latest
homo_sapiens_core.meta_coord file.  Then it will update the meta_coord
table for all the following table names one by one
  regulatory_feature
  probe_feature
  external_feature
  annotated_feature
 
  );

if (scalar @ARGV == 0 ) {
  print $usage, "\n";
  exit 0;
}

GetOptions( 'help'        => \$help,
            'host|h=s'      => \$host,
            'port=i'      => \$port,
            'user|u=s'      => \$user,
            'pass|p=s'      => \$pass,
            'dbpattern|d=s' => \$dbpattern );

if ($help) {
  print $usage, "\n";
  exit 0;
}

#print "help: $help argv:"
#  . scalar(@ARGV)
#  . "$host $port $user $pass $dbpattern\n";

my $dsn = "DBI:mysql:host=$host";
$dsn .= ";port=$port" if ($port);

my $db = DBI->connect( $dsn, $user, $pass ) || die 'Failed to connect to DB';


my @dbnames =
  map { $_->[0] } @{ $db->selectall_arrayref("show databases") };

my $found = 0;

for my $dbname (@dbnames) {

  next if ( $dbname !~ /$dbpattern/ );
  $found = 1;

  print "Updating $dbname meta_coord\n";

  my $dbc =
    new Bio::EnsEMBL::DBSQL::DBConnection( -host   => $host,
                                           -port   => $port,
                                           -user   => $user,
                                           -pass   => $pass,
                                           -dbname => $dbname );

  my @table_names = qw(
					   regulatory_feature
					   probe_feature
					   external_feature
					   annotated_feature
					  );

  unless (
           system(     "mysql -h$host -P$port -u$user -p$pass -N "
                     . "-e 'SELECT * FROM meta_coord' $dbname "
                     . "> $dbname.meta_coord.backup"
           ) == 0 )
  {
    print "Can't dump the original meta_coord for back up\n";
    exit 1;
  } else {
    print "Original meta_coord table backed up in "
      . "$dbname.meta_coord.backup\n";
  }

  foreach my $table_name (@table_names) {
    print "Updating $table_name table entries...";
    my $sql = "DELETE FROM meta_coord WHERE table_name = ?";
    my $sth = $dbc->prepare($sql);
    $sth->execute($table_name);
    $sth->finish;

    $sql =
        "INSERT INTO meta_coord "
      . "SELECT '$table_name', s.coord_system_id, "
      . "MAX( t.seq_region_end - t.seq_region_start + 1 ) "
      . "FROM $table_name t, seq_region s "
      . "WHERE t.seq_region_id = s.seq_region_id "
      . "GROUP BY s.coord_system_id";
    $sth = $dbc->prepare($sql);
    $sth->execute;
    $sth->finish;
    print "Done\n";
  }
} ## end for my $dbname (@dbnames)

die "Found no databases matching pattern /$dbpattern/\n" if ! $found;
