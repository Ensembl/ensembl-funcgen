#!/usr/local/ensembl/bin/perl -w


use strict;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::HealthChecker;

my $reg = "Bio::EnsEMBL::Registry";

my $species = 'mus_musculus';
my $schema_build = '48_37a';
my $pass = shift @ARGV;
my $port = 3306;
my $user = 'ensadmin';
my $host = 'ens-genomics1';
my $skip_mc;
my @builds =('');#blank string for the default build;
push @builds, @ARGV;

#only loads v43 no v44 DBs???
#$reg->load_registry_from_db(
#							-host => 'ens-staging',
#							-user => 'ensro',
#							#-port => '3362',
#							-verbose => 1,
#						   );#







my $dnadb = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
												-dbname  => "${species}_core_${schema_build}",
												-host    => 'ens-staging',# 
												#-host    => 'ensembldb',
												-port    => 3306,
												-user    =>  $user,
												#-user    =>  'anonymous',
												-pass    => $pass,
												-species => $species,
											   );
die ("You have not provided sufficient arguments to make a DB connection\n".
	 "-dbname  => ${species}_funcgen_${schema_build}, host   => $host, port   => $port, user   => $user, pass   => $pass") if ! $dnadb;


#This will add the default chromosome CS, but not any other levels
my $efg_db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
														  -dbname  => "${species}_funcgen_${schema_build}",
														  -host    => $host,
														  -port    => $port,
														  -user    => $user,
														  -pass    => $pass,
														  -species => $species,
														  -dnadb   => $dnadb,
														 );


my $hchecker = Bio::EnsEMBL::Funcgen::Utils::HealthChecker->new(
																-db     => $efg_db,
																-builds => \@builds,
															   );

$hchecker->update_db_for_release($skip_mc);
