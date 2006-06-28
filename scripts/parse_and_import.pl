#!/opt/local/bin/perl -w


=head1 NAME

ensembl-efg parse_and_import.pl
  
=head1 SYNOPSIS

parse_and_import.pl [options]

Options:

Mandatory
  -instance|i      Instance name
  -format|f        Data format
  -group|g         Group name

Optional
  -pass|p          The password for the target DB, if not defined in GroupDefs.pm
  -data_root       The root data dir
  -fasta           Fasta dump flag
  -norm|n          Normalisation method (default=vsn)
  -species|s       Species name(latin e.g. "Homo sapiens" or "Mus musculus")
  -debug|d         Debug level (1-3)
  -log_file|l      Defines the log file
  -help            Brief help message
  -man             Full documentation

=head1 OPTIONS

=over 8

=item B<-instance|i>

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

=item B<-norm|n>

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
use Bio::EnsEMBL::Funcgen::Experiment;

use strict;



my ($input_name, $instance, $output_dir, $group, $pass, $help, $man, $species, $nmethod);


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
			"instance|i=s" => \$instance,
			"format|f=s"   => \$format,
			"vendor|v=s"   => \$vendor,
			"pass|p=s"     => \$pass,
			"group|g=s"    => \$group,#Need user here too? Use group defs to avoid typos?
			"species|s=s"  => \$species,
			"debug|d=i"    => \$main::_debug_level,
			"data_root=s"  => \$data_dir,
			"fasta"        => \$fasta,
			"recover|r"    => \$recover,
			"norm|n=s"     => \$nmethod,
			"log_file|l=s" => \$main::_log_file,
			"debug_file=s" => \$main::_debug_file,
			#should have MAGE flag here? or would this be format?
			"interactive"  => \$interactive,
			"help|?"       => \$help,
			"man|m"        => \$man,			
		   );

pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

#Build and validate all these in Experiment::new? We only need these for importing/analysing???
$output_dir  = $data_dir."/".$vendor."/".$instance;
$main::_log_file = $output_dir."/${instance}.log" if(! defined $main::_log_file);
$main::_debug_file = $output_dir."/${instance}.dbg" if(! defined $main::_debug_file);


### CREATE EXPERIMENT AND VALIDATE ###
#Write Importer instead?
#And trim down Experiment?

my $Exp = Bio::EnsEMBL::Funcgen::Experiment->new(
												 instance    => $instance,
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
												 #Exp does not build input dir, but could
												 #This allows input dir to be somewhere 
												 #other than efg dir structure
												);


#Validate, parse and import all experiment data

$Exp->register_experiment();


1;
