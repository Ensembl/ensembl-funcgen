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
my ($exp_name, $help, $pids, $man, @features, $chr, $not_status);

my $dbhost = $ENV{'EFG_HOST'};
my $port = $ENV{'EFG_PORT'};
my $user = $ENV{'EFG_READ_USER'};
my $dbname = $ENV{'EFG_DBNAME'};
my $cdbname = $ENV{'CORE_DBNAME'};
my $no_zip = 0;

my $anal_name = 'Nessie';
my $out_dir = ".";

GetOptions (
			"feature_set=s"    => \$fset_name,
			"pass=s"           => \$pass,
			"port=s"           => \$port,
            "dbname=s"         => \$dbname,
			"dbhost=s"         => \$dbhost,
            "cdbname=s"        => \$cdbname,
			"outdir=s"         => \$out_dir,
			'no_zip|z'         => \$no_zip,
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

throw("Must pass a feature_set name via -feature_set, i.e. 'RegulatoryFeatures'") if ! $fset_name;

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
my ($outline, @output);
my $fset_ftype = ucfirst($fset->type).'Feature';

foreach my $slice(@{$slice_a->fetch_all('toplevel')}){
  my $cnt = 0;

  print "Dumping slice ".$slice->name."\n";

  my $chr_name = 'Chr'.$slice->seq_region_name();

  my $ofile = $out_dir."/".$fset_name.'.'.$chr_name.'.gff';
  
  open (OUT, ">$ofile") || die ("Can't open $ofile for writing");

  foreach my $feature(@{$fset->get_Features_by_Slice($slice)}){
	$cnt++;
	
	#seqid source type start end score strand phase attrs

	#Let's push this onto an array here and only print OUT every 5000 lines to save I/O?
	#We are no handling strand here

	$outline = join("\t", ($chr_name, $dbname, $fset_ftype, $feature->start(), $feature->end(), '.', '.', '.', 'Name='.$feature->feature_type->name.';'));

	#associated_feature_types?
	#http://www.sequenceontology.org/gff3.shtml

	if($fset_ftype eq 'RegulatoryFeature'){
	  $outline .= join('; ', (' ID='.$feature->stable_id(), 'Note=Consists of following features: '.join(',', map {join(':', $_->feature_type->name, $_->cell_type->name)}@{$feature->regulatory_attributes()})));

	}
	#elsif($fset_ftype eq 'AnnotatedFeature'){
	#  #feature_type->name
	#  
	#}				
	elsif($fset_ftype eq 'ExternalFeature'){
	  #It may be more appropriate to have this as the Name for external_features?
	  $outline .= ' Alias='.$feature->display_label;
	}

	push @output, $outline;


	if(scalar(@output) == 1000){
	  print OUT join("\n", @output)."\n";
	  @output = ();
	}
  }

  print OUT join("\n", @output);
  @output = ();

  close(OUT);

  #Remove empty files and zip

  if(-z $ofile){
	print "No features found\n";
	unlink $ofile;
  }
  else{
	print "Found $cnt features\n";

	if(! $no_zip){
	  system("gzip $ofile");
	}
  }
}
