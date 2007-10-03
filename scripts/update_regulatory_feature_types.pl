#!/usr/local/ensembl/bin/perl -w


use strict;
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils;

my $reg = "Bio::EnsEMBL::Registry";

my $species = 'homo_sapiens';
my $schema_build = '47_36i';
my $pass = shift @ARGV;
my $port = 3306;
my $user = 'ensadmin';
my $host = 'ens-genomics1';
my $out_dir = './';
my $log_file = $out_dir.'RegulatoryFeatures.type.log';

my %regex_types =();


my $efg_db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
														  -dbname  => $species.'_funcgen_'.$schema_build,
														  -host    => $host,
														  -port    => $port,
														  -user    => $user,
														  -pass    => $pass,
														  -species => $species,
														 );


die ("You have not provided sufficient arguments to make a DB connection\n".
	 "-dbname  => ${species}_funcgen_{$schema_build}, host   => $host, port   => $port, user   => $user, pass   => $pass") if ! $efg_db;


my $ftype_adaptor = $efg_db->get_FeatureTypeAdaptor();
my %ftypes = ();

foreach my $type(values %regex_types){

  my $ftype = $ftype_adaptor->fetch_by_name($type);

  throw('FeatureType $type not present in the DB, use import_type.pl before running update') if ! defined $ftype;

  $ftypes{$type} = $ftype;
}

my $handle = open_file($log_file, '>');

my $rf_fset = $efg_db->get_RegulatoryFeatureAdaptor()->fetch_by_name('RegulatoryFeatures');

#should we split this up into slices?
my $sql1 = 'UPDATE regulatory_feature where regulatory_feature_id=';
my $sql2 = ' set feature_type_id=';
my $cnt = 0;

my @features = @{$rf_fset->fetch_all_by_FeatureSet($rf_fset)};

print 'Assigning new feature types: '.join("\t", values(%regex_types))."\n";
print 'Assigning to '.scalar(@features)." RegulatoryFeatures\n";

FEATURE: foreach my $rfeat(@features){
  my $assigned = 0;

  foreach my $regex(keys %regex_types){


	#access directly just in case we put a method hack inplace to retain the binary string, but display something moire meaningful
	if($rfeat->{'display_label'} =~ /$regex/){

	  if($assigned){
		print $handle 'WARNING: Skipping RegulatoryFeature with multiple type assignments:\tdbID:'.$rfeat->dbID."\t".$rfeat->feature_type->name().' & '.$regex_types{$regex}."\n";
		next FEATURE;
	  }
	  
	  $rfeat->feature_type($ftypes{$regex_types{$regex}});
	}
  }

  if(! $assigned){
	print $handle "WARNING: Skipping unassigned RegulatoryFeature:\tdbID:".$rfeat->dbID()."\n";
  }
  else{
	$cnt ++;
	$efg_db->dbc->db_handle->do($sql1.$rfeat->dbID().$sql2.$rfeat->feature_type->dbID());
  }
  
  print "Updated $cnt RegulatoryFeatures\n" if(! ($cnt  % 10000));

}
print "Finished updating a total of $cnt / ".scalar(@features)." RegulatoryFeature\n";
