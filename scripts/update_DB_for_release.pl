
use strict;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::HealthChecker;
use Pod::Usage;
use Getopt::Long;

my ($pass, $species, $schema_build, $skip_meta_coord, $dnadb_host, $dnadb, $check_displayable);
my ($help, $man, $dbname);
my $user = 'ensadmin';
my $port = 3306;
my $host = 'ens-genomics1';
my @tmp_args = @ARGV;

GetOptions( 
		   'host|h=s'         => \$host,
		   'port=i'           => \$port,
		   'user|u=s'         => \$user,
		   'pass|p=s'         => \$pass,
		   'dbname=s'         => \$dbname,
		   'species|d=s'      => \$species,
		   'data_versions=s' => \$schema_build,#mandatory as default ensembldb will not exist
		   'dnadb_host=s'     => \$dnadb_host,
		   'skip_meta_coord'  => \$skip_meta_coord,
		   'check_displayable' => \$check_displayable,
		   'help|?'             => \$help,
		   'man|m'              => \$man,
		   'log_file=s'        => \$main::_log_file,
		   'tee'               => \$main::_tee,
		   #'slice=s'          => \$test_slice,
		   #skip dumps?
		   #force update
		   #add opt for old, new & stable fset name
		  ) or pod2usage(
						 -exitval => 1,
						 -message => "Params are:\t@tmp_args"
						);

pod2usage(0) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;


my @builds = @ARGV;

$dbname = "${species}_funcgen_${schema_build}" if ! defined $dbname;
$main::_log_file ||= $ENV{'HOME'}."/logs/update_DB_for_release.$dbname.$$.log";
print "Writing log to:\t".$main::_log_file."\n";

if($dnadb_host){
  $dnadb = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
											   -dbname  => "${species}_core_${schema_build}",
											   -host    => $dnadb_host,
											   #-host => 'ensdb-1-13',
											   #-port => 5307,
											   -port    => 3306,
											   -user    =>  'ensro',
											   #-pass    => $pass,
											   -species => $species,
											   -group   => 'core',
											   );
}


#This will add the default chromosome CS, but not any other levels
my $efg_db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
														  -dbname  => $dbname,
														  -host    => $host,
														  -port    => $port,
														  -user    => $user,
														  -pass    => $pass,
														  -species => $species,
														  -dnadb   => $dnadb,
														 );



#Test the db connections
$efg_db->dbc->db_handle;
$efg_db->dnadb->dbc->db_handle;

my $hchecker = Bio::EnsEMBL::Funcgen::Utils::HealthChecker->new(
																-db     => $efg_db,
																-builds => \@builds,
																-skip_meta_coord => $skip_meta_coord,
																#-force_update
																#skip_dumps
																-check_displayable => $check_displayable,
															   );

$hchecker->update_db_for_release();
$hchecker->log_data_sets();
$hchecker->check_stable_ids;
