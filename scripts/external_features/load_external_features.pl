# Load data from a file of regulatory regions into a database

#To do
#remove get_all_regulatory_features from Gene.pm
#do we need some different feature tables for these high volume features?
#how would we trace this feature_type association...in feature_set, feature_type column?
#would be feature_type for most, maybe rna_feature and/or motif_feature
#store influence and evidence in a functional xref table?
#set analysis and ftype adaptor internally?
#implement version in xref to denote version of stable id and hence implicitly release of core DB
#Or implement release in external_db?
#can we use linkage annotation in object_xref?

use strict;

use Bio::EnsEMBL::Funcgen::Parsers::vista;
use Bio::EnsEMBL::Funcgen::Parsers::cisred;
use Bio::EnsEMBL::Funcgen::Parsers::miranda;
use Bio::EnsEMBL::Funcgen::Parsers::redfly;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Getopt::Long;

my ($host, $user, $pass, $port, $dbname, $species, $file, $clobber, $type, $old_assembly, $new_assembly, $archive);

GetOptions( "host=s",         \$host,
			"user=s",         \$user,
			"species=s",      \$species,
			"pass=s",         \$pass,
			"port=i",         \$port,
			"dbname=s",       \$dbname,
			"type=s",         \$type,
			"file=s",         \$file,
			"archive=s",      \$archive,
			"clobber",        \$clobber,#This is not behaving like a boolean flag??
			"help",           \&usage,
			"old_assembly=s", \$old_assembly,
			"new_assembly=s", \$new_assembly
	  );

usage() if (!$host || !$user || !$dbname || !$type);




my $cdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
											  -host   => 'ens-staging',
											  -port   => $port,
											  -user   => $user,
											  -pass   => $pass,
											  -species => $species,
											  -group  => 'core',
											  -dbname => $species.'_core_49_44',
											 );


my $db_adaptor = new Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor(
															 -host   => $host,
															 -port   => $port,
															 -user   => $user,
															 -pass   => $pass,
															 -species => $species,
															 -group  => 'funcgen',
															 -dbname => $dbname,
															 -dnadb  => $cdb,
															);

#test db connections
$db_adaptor->dbc->db_handle;
$db_adaptor->dnadb->dbc->db_handle;
$type = lc($type);
#we need to validate module in a different way
# validate type
#exit(1) if (!RegulatoryFeatureParser::BaseParser::validate_type($db_adaptor, $type));
#Need to eval the use? No use is performed before anything?


eval "require Bio::EnsEMBL::Funcgen::Parsers::$type";
if($@) {
  die("Did not find a parser module corresponding to $type");
}

my $parser = "Bio::EnsEMBL::Funcgen::Parsers::$type"->new(
												   -db      => $db_adaptor,
												   -clobber => $clobber,
												   -archive => $archive,
												   -type    => $type,
												  );


#RegulatoryFeatureParser::BaseParser::delete_existing($db_adaptor, $type) if ($del);

#if(! defined $file){
#  print "Must provide a file to parse\n";
#  &usage;
#}

#if(! -e $file){
 # print "Can't find $file\n";
 # &usage;
#}

$parser->parse_and_load($file, $old_assembly, $new_assembly);

#$parser->upload_features_and_factors($objects);




sub usage {
  print <<EOF; exit(0);

Usage: perl $0 <options>

  -host    Database host to connect to

  -port    Database port to connect to

  -dbname  Database name to use

  -user    Database username

  -pass    Password for user

  -type    Type of regulatory features to upload; must be present in analysis table and
           have a correspoinding parser e.g. vista, cisred, miranda.

  -file    File to load data from; must be parseable by parser corresponding to type

  -clobber Delete all features from type feature_set(s) before loading

  -archive Specify a version suffix to move the old feature_set to e.g. 47 will move feature_set to feature_set_name_v47

  -old_assembly, -new_assembly If *both* these options are specified, co-ordinates will be
           projected from the old assembly to the new assembly before being stored. This
           relies on a suitable mapping between the two assemblies having already been 
           carried out.

  -help    This message

EOF

}

