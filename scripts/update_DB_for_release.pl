#!/usr/local/ensembl/bin/perl -w

use strict;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::HealthChecker;
use Pod::Usage;
use Getopt::Long;

my ($pass, $species, $schema_build, $skip_meta_coord, $from_ensdb, $dnadb, $check_displayable);
my ($help, $man, $dbname);
my $user = 'ensadmin';
my $port = 3306;
my $host = 'ens-genomics1';

GetOptions( 'host|h=s'         => \$host,
            'port=i'           => \$port,
            'user|u=s'         => \$user,
            'pass|p=s'         => \$pass,
			'dbname=s'         => \$dbname,
			#'slice=s'          => \$test_slice,
            'species|d=s'      => \$species,
			'data_versions=s' => \$schema_build,#mandatory as default ensembldb will not exist
			'from_ensembldb'   => \$from_ensdb,
			'skip_meta_coord'  => \$skip_meta_coord,
			#skip dumps?
			#force update
			'check_displayable' => \$check_displayable,
			"help|?"             => \$help,
			"man|m"              => \$man,
			#add opt for old, new & stable fset name
		   );


pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;



my @builds = @ARGV;

$dbname = "${species}_funcgen_${schema_build}" if ! defined $dbname;

if(! $from_ensdb){
  $dnadb = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
											   -dbname  => "${species}_core_${schema_build}",
											   -host    => 'ens-staging',
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

