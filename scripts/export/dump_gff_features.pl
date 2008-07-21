#!/usr/local/ensembl/bin/perl

use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::ResultSet;
use Bio::EnsEMBL::Utils::Exception qw( throw );

$| =1;

my ($file, $ofile, $pass, $line, $fset_name);
my ($exp_name, $dbname, $cdbname, $help, $pids, $man, @features, $chr, $not_status);
my $dbhost = $ENV{'EFG_HOST'};
my $port = $ENV{'EFG_PORT'};
my $user = $ENV{'EFG_READ_USER'};
my $anal_name = 'Nessie';
my $out_dir = ".";


GetOptions (
			"feature_set=s" => \$fset_name,
			"pass=s"           => \$pass,
			"port=s"           => \$port,
            "dbname=s"         => \$dbname,
			"dbhost=s"         => \$dbhost,
            "cdbname=s"        => \$cdbname,
			"outdir=s"         => \$out_dir,
			"help|?"           => \$help,
			"man|m"            => \$man,
		   );

pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;


### Set up adaptors and FeatureSet 

if(! $cdbname  || ! $dbname ){
  throw("You must provide a funcgen(-dbname) and a core(-cdbname) dbname");
}
throw("Must define your funcgen dbhost -dbhost") if ! $dbhost;
#throw("Must supply an input file with -file") if ! $file;
#throw("Must supply a password for your fungen db") if ! $pass;

#this need genericising for ensembldb/ens-livemirror

my $cdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
											  -host => "ens-staging",
											  -dbname => $cdbname,
											  #-species => "homo_sapiens",
											  -user => $user,
											  -pass => "",
											  -group => 'core',
											  -port => 3306,
											 );

my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
													  -host => $dbhost,
													  -dbname => $dbname,
													  #-species => "homo_sapiens",
													  -user => $user,
													  -pass => $pass,
													  -dnadb => $cdb,
													  -port => $port,
													 );


#Check DB connections
$cdb->dbc->db_handle;
$db->dbc->db_handle;

#Check out_dir
#Use Helper!
system('mkdir -p '.$out_dir) if(! -d $out_dir);


#ChrX  . operon   XXXX YYYY  .  +  . ID=operon01;name=my_operon

my $slice_a = $db->get_SliceAdaptor();
my $fset_a = $db->get_FeatureSetAdaptor();
my $fset = $fset_a->fetch_by_name($fset_name);

foreach my $slice(@{$slice_a->fetch_all('toplevel')}){

  print "Dumping slice $slice ".$slice->seq_region_name."\n";

  my $chr_name = 'Chr'.$slice->seq_region_name();

  my $ofile = $out_dir."/".$fset_name.'.'.$chr_name.'.gff';
  
  open (OUT, ">$ofile") || die ("Can't open $ofile for writing");

  foreach my $feature(@{$fset->get_Features_by_Slice($slice)}){
	
	#seqid source type start end score strand phase attrs

	print OUT join("\t", ('>'.$chr_name, $dbname, 'regulatory feature', 
                      $feature->start(), $feature->end(), '.', '.', '.', 'ID='.$feature->stable_id().';Note=Consists of following features: '.
                      join(',', map {join(':', $_->feature_type->name, $_->cell_type->name)}@{$feature->regulatory_attributes()})))."\n";

  }

  close(OUT);

}
