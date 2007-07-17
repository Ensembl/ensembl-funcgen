#!/usr/local/ensembl/bin/perl -w


use strict;
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::ProbeFeature;

my $reg = "Bio::EnsEMBL::Registry";

my $species = 'homo_sapiens';
my $schema_build = '46_36h';
my $pass = $ARGV[0];
my $port = 3306;
my $user = 'ensadmin';
my $host = 'ens-genomics1';

#only loads v43 no v44 DBs???
#$reg->load_registry_from_db(
#							-host => 'ens-staging',
#							-user => 'ensro',
#							#-port => '3362',
#							-verbose => 1,
#						   );#







my $dnadb = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
												-dbname  => $species.'_core_'.$schema_build,
												-host    => 'ens-staging',
												-port    => 3306,
												-user    => $user,
												-pass    => $pass,
												-species => $species,
											   );
die ("You have not provided sufficient arguments to make a DB connection\n".
	 "-dbname  => ${species}_funcgen_{$schema_build}, host   => $host, port   => $port, user   => $user, pass   => $pass") if ! $dnadb;


#This will add the default chromosome CS, but not any other levels
my $efg_db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
														  -dbname  => $species.'_funcgen_'.$schema_build,
														  -host    => $host,
														  -port    => $port,
														  -user    => $user,
														  -pass    => $pass,
														  -species => $species,
														  -dnadb   => $dnadb,
														 );


die ("You have not provided sufficient arguments to make a DB connection\n".
	 "-dbname  => ${species}_funcgen_{$schema_build}, host   => $host, port   => $port, user   => $user, pass   => $pass") if ! $efg_db;



my $core_db = $efg_db->dnadb();


#do we need to add the none default levels here?
#or are we only bothered about those which constitute the toplevel?

#To make sure we have all the correct levels in eFG we need to get all the names.
#then get all by name from the core db and set them as the dnadb.
# we also need to get all the toplevel seq_regions and store them in the seq_region table
#use BaseFeatureAdaptor::_pre_store with and array of pseudo feature on each top level slice

my $pf_adaptor = $efg_db->get_ProbeFeatureAdaptor();
my $slice_adaptor = $efg_db->get_SliceAdaptor();

foreach my $slice(@{$slice_adaptor->fetch_all('toplevel')}){

  if($slice->start() != 1){
	#we must have some sort of PAR linked region i.e. Y
	$slice = $slice_adaptor->fetch_by_region($slice->coord_system_name(), $slice->seq_region_name());
  }

  print "Storing seq_region info for slice:\t".$slice->name()."\n";

  my $pseudo_feature = Bio::EnsEMBL::Funcgen::ProbeFeature->new
	(
	 -slice => $slice,
	 -start => $slice->start(),
	 -end   => $slice->end(),
	 -strand => 0,
	);

  $pf_adaptor->_pre_store($pseudo_feature);
}
