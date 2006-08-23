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
  -species|s       Species name any standard e! species alias(will be reset to dbname/latin name e.g. "homo_sapiens")
  -location        Physical location of experimental group
  -contact         Contact details for experimental group
  -debug           Debug level (1-3)
  -log_file        Defines the log file
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
use Bio::EnsEMBL::Funcgen::Importer;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );

use strict;



my ($input_name, $name, $output_dir, $loc, $contact, $group, $pass);
my ($data_version, $help, $man, $species, $nmethod, $dnadb);
my $reg = "Bio::EnsEMBL::Registry";

#to be removed
my ($input_dir, $import_dir);

my $data_dir = $ENV{'EFG_DATA'};
my $interactive = 1;
my $format = "tiled";
my $vendor = "nimblegen";
my $user = "ensadmin";
my $host = 'localhost';
my $port = '3306';
my $fasta = 0;#Shouldn't this be on by default?
my $recover = 0;
my $verbose = 0;

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
			"port|l=s"     => \$port,
			"host|h=s"     => \$host,
			"user|u=s"     => \$user,
			"group|g=s"    => \$group,#Need user here too? Use group defs to avoid typos?
			"species|s=s"  => \$species,
			"data_version|d=s" => \$data_version,
			"debug=i"    => \$main::_debug_level,
			"data_root=s"  => \$data_dir,
			"fasta"        => \$fasta,
			"recover|r"    => \$recover,
			"norm=s"       => \$nmethod,
			"location=s"   => \$loc,
			"contact=s"    => \$contact,
			"log_file=s"   => \$main::_log_file,
			"debug_file=s" => \$main::_debug_file,
			#should have MAGE flag here? or would this be format?
			"interactive"  => \$interactive,
			"help|?"       => \$help,
			"man|m"        => \$man,
			"verbose=s"    => \$verbose,
		   );


#Need to add (primary) design_type and description, or add to defs file?


pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

#Build and validate all these in Experiment::new? We only need these for importing/analysing???
$output_dir  = $data_dir."/".$vendor."/".$name;
mkdir $output_dir;#log/debug files fail in Helper without this
$main::_log_file = $output_dir."/${name}.log" if(! defined $main::_log_file);
$main::_debug_file = $output_dir."/${name}.dbg" if(! defined $main::_debug_file);



### SET UP IMPORTER (FUNCGENDB/DNADB/EXPERIMENT) ###

my $Imp = Bio::EnsEMBL::Funcgen::Importer->new(
											   name        => $name,
											   format      => $format,
											   vendor      => $vendor,
											   group       => $group,
											   pass        => $pass,
											   host        => $host,
											   user        => $user,
											   port        => $port,

											   
											   data_version => $data_version,
											   data_root   => $data_dir,
											   output_dir  => $output_dir,
											   recover     => $recover,
											   dump_fasta  => $fasta,
											   norm_method => $nmethod,
											   species     => $species,
											   location    => $loc,
											   contact     => $contact,
											   verbose     => $verbose,
											   #Exp does not build input dir, but could
											   #This allows input dir to be somewhere 
											   #other than efg dir structure
											  );



#Can be moved to set up script?

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


my $t_anal = Bio::EnsEMBL::Analysis->new(
									   -logic_name      => 'TilingHMM',
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
									   -description    => 'HMM based predictions based on tiling data',
									   -display_label   => 'TilingHMM',
									  );


#this checks if already stored
$anal_a->store($anal);
$anal_a->store($raw_anal);
$anal_a->store($vsn_anal);
$anal_a->store($t_anal);


#Validate, parse and import all experiment data


$Imp->register_experiment();


1;
