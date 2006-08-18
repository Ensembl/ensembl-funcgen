#!/opt/local/bin/perl -w


=head1 NAME

ensembl-efg parse_and_import.pl
  
=head1 SYNOPSIS

parse_and_import.pl [options]

Options:

Mandatory
  -name|n          Instance name
  -format|f        Data format
  -group|g         Group name

Optional
  -pass|p          The password for the target DB, if not defined in GroupDefs.pm
  -data_root       The root data dir
  -fasta           Fasta dump flag
  -norm|n          Normalisation method (default=vsn)
  -species|s       Species name(latin e.g. "Homo sapiens" or "Mus musculus")
  -location        Physical location of experimental group
  -contact         Contact details for experimental group
  -debug|d         Debug level (1-3)
  -log_file|l      Defines the log file
  -help            Brief help message
  -man             Full documentation

=head1 OPTIONS

=over 8

=item B<-name|n>

Mandatory:  Instance name for the data set, this is the directory where the native data files are located

=item B<-format|f>

Mandatory:  The format of the data files e.g. nimblegen

=over 8

=item B<-group|g>

Mandatory:  The name of the experimental group

=over 8

=item B<-data_root>

The root data dir containing native data and pipeline data, default = $ENV{'EFG_DATA'}

=over 8

=item B<-fasta>

Flag to turn on dumping of all probe_features in fasta format for the remapping pipeline

=item B<-norm>

Normalisation method, deafult is the Bioconductor vsn package which performs generalised log ratio transformations

=item B<-species|s>

Species name for the array.

=item B<-debug>

Turns on and defines the verbosity of debugging output, 1-3, default = 0 = off

=over 8

=item B<-log_file|l>

Defines the log file, default = "${instance}.log"

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> take several options, including an definitions file to parse and import array data into the ensembl-efg DB

=cut


#add @INC stuff here, or leave to .bashrc/.efg?

BEGIN{
	if(! defined $ENV{'EFG_DATA'}){
		if(-f "~/src/ensembl-efg/scripts/.efg"){
			system (". ~/src/ensembl-efg/scripts/.efg");
		}else{
			die ("This script requires the .efg file available from ensembl-efg\n".
				 "Please source it before running this script\n");
		}
	}
}
	

#use Bio::EnsEMBL::Root; #Only used for rearrange see pdocs
#Roll own Root object to handle debug levels, logging, dumps etc.

### MODULES ###
use Getopt::Long;
#use Carp;#For dev only? cluck not exported by default Remove this and implement in Helper
use Pod::Usage;
#POSIX? File stuff
#use Data::Dumper;
use Bio::EnsEMBL::Funcgen::Importer;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Registry;


use strict;



my ($input_name, $name, $output_dir, $loc, $contact, $group, $pass, $help, $man, $species, $nmethod);
my $reg = " Bio::EnsEMBL::Registry";

#to be removed
my ($input_dir, $import_dir);

my $data_dir = $ENV{'EFG_DATA'};
my $interactive = 1;
my $format = "tiled";
my $vendor = "nimblegen";
my $fasta = 0;#Shouldn't this be on by default?
my $recover = 0;

#Definitely need some sort of Defs modules for each array?

$main::_debug_level = 0;

my @steps = ("s1_init_exp");#do we need steps?

#Use some sort of DBDefs for now, but need  to integrate with Register, and have put SQL into (E)FGAdaptor?
#Use ArrayDefs.pm module for some of these, class, vendor, format?
#ArrayDefs would also contain paths to data and vendor specific parse methods?

GetOptions (
			"name|n=s"     => \$name,
			"format|f=s"   => \$format,
			"vendor|v=s"   => \$vendor,
			"pass|p=s"     => \$pass,
			"group|g=s"    => \$group,#Need user here too? Use group defs to avoid typos?
			"species|s=s"  => \$species,
			"debug|d=i"    => \$main::_debug_level,
			"data_root=s"  => \$data_dir,
			"fasta"        => \$fasta,
			"recover|r"    => \$recover,
			"norm=s"       => \$nmethod,
			"location=s"   => \$loc,
			"contact=s"    => \$contact,
			"log_file|l=s" => \$main::_log_file,
			"debug_file=s" => \$main::_debug_file,
			#should have MAGE flag here? or would this be format?
			"interactive"  => \$interactive,
			"help|?"       => \$help,
			"man|m"        => \$man,			
		   );


#Need to add (primary) design_type and description, or add to defs file?


pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

#Build and validate all these in Experiment::new? We only need these for importing/analysing???
$output_dir  = $data_dir."/".$vendor."/".$name;
mkdir $output_dir;#log/debug files fail in Helper without this
$main::_log_file = $output_dir."/${name}.log" if(! defined $main::_log_file);
$main::_debug_file = $output_dir."/${name}.dbg" if(! defined $main::_debug_file);


### CREATE EXPERIMENT AND VALIDATE ###
#Write Importer instead?
#And trim down Experiment?

#Need to generate and pass dbadaptor here, else use default settings(harcoded at mo)
#retrieve v25 DB for current dataset


$reg->load_registry_from_db(-host => 'localhost',
							-user => 'ensdmin',
							-pass => 'ensembl',
						   );


my $db = $reg->get_DBAdaptor('Human', 'funcgen');


print "Got DBAdaptor $db :\n".Data::Dumper::Dumper(\$db)."\n";

exit;

							

my $Imp = Bio::EnsEMBL::Funcgen::Importer->new(
											   name        => $name,
											   format      => $format,
											   vendor      => $vendor,
											   group       => $group,
											   pass        => $pass,
											   data_root   => $data_dir,
											   output_dir  => $output_dir,
											   recover     => $recover,
											   dump_fasta  => $fasta,
											   norm_method => $nmethod,
											   species     => $species,
											   location    => $loc,
											   contact     => $contact,
											   #Exp does not build input dir, but could
											   #This allows input dir to be somewhere 
											   #other than efg dir structure
												);


#so we need to specify a dnadb for the register epxeriment step

my $anal_a = $Imp->db->get_AnalysisAdaptor();
my $anal = Bio::EnsEMBL::Analysis->new(
									   -logic_name      => 'VendorMap',
									   -db              => 'NULL',
									   -db_version      => 'NULL',
									   -db_file         => 'NULL',
									   -program         => 'NULL',
									   -program_version => 'NULL',
									   -program_file    => 'NULL',
									   -gff_source      => 'NULL',
									   -gff_feature     => 'NULL',
									   -module          => 'NULL',
									   -module_version  => 'NULL',
									   -parameters      => 'NULL',
									   -created         => 'NULL',
									   -description    => 'Original mapping provided by Array Vendor',
									   -display_label   => 'Vendor mapping',
									  );


my $raw_anal = Bio::EnsEMBL::Analysis->new(
										   -logic_name      => 'RawValue',
										   -db              => 'NULL',
										   -db_version      => 'NULL',
										   -db_file         => 'NULL',
										   -program         => 'NULL',
										   -program_version => 'NULL',
										   -program_file    => 'NULL',
										   -gff_source      => 'NULL',
										   -gff_feature     => 'NULL',
										   -module          => 'NULL',
										   -module_version  => 'NULL',
										   -parameters      => 'NULL',
										   -created         => 'NULL',
										   -description    => 'Raw value',
										   -display_label   => 'Raw value',
										  );

my $vsn_anal = Bio::EnsEMBL::Analysis->new(
										   -logic_name      => 'VSN',
										   -db              => 'NULL',
										   -db_version      => 'NULL',
										   -db_file         => 'NULL',
										   -program         => 'NULL',
										   -program_version => 'NULL',
										   -program_file    => 'NULL',
										   -gff_source      => 'NULL',
										   -gff_feature     => 'NULL',
										   -module          => 'NULL',
										   -module_version  => 'NULL',
										   -parameters      => 'NULL',
										   -created         => 'NULL',
										   -description    => 'Generalised log transformation',
										   -display_label   => 'VSN',
										  );


#print "anal adaptor is $anal_a dump:\n".Data::Dumper::Dumper(\$anal_a).'\n';


#this checks if already stored
$anal_a->store($anal);
$anal_a->store($raw_anal);
$anal_a->store($vsn_anal);


#Validate, parse and import all experiment data

#hardcoded..also hardcoded in DBAdaptor autogeneration, need to validation species name, store in meta.species_latin?
#my $species = "homo_sapiens";
#my $schema_build = "25_34e";#This is the true schema_build for the first Nimblegen set.

#my $dnadb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
#		  (
#		   -host => "ensembldb.ensembl.org",
#		   -user => "anonymous",
#		   -dbname => "${species}_core_${schema_build}",
#		  );

#Need to generate cood_system from dnadb, should need this if we're importing vie store i.e. on a slice.
#But will for flat file

	

$Imp->register_experiment();


1;
