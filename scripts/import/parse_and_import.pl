#!/software/bin/perl -w


=head1 NAME

 ensembl-efg parse_and_import.pl
  
=head1 SYNOPSIS

 parse_and_import.pl [options]

=head1 OPTIONS

 Options:

 Mandatory
  -name|n          Instance name
  -format|f        Data format
  -group|g         Group name

 Optional
  -registry_host Host to load registry from(defaults to ensembldb.ensembl.org)
  -registry_port Port for registry host
  -registry_user User for registry host 
  -registry_pass Password for registry user
  -pass|p          The password for the target DB, if not defined in GroupDefs.pm
  -data_root       The root data dir
  -dbname          Defines the eFG dbname if it is not standard
  -ssh             Forces use of tcp protocol for using ssh forwarded ports from remote machine(remove and just specify 127.0.0.1 as host?)
  -array_set       Flag to treat all chip designs as part of same array
  -array_name      Name of the set of array chips
  -result_set      Name to give the raw/normalised result set.
  -result_files    Space separated list of result files paths
  -fasta           Fasta dump flag
  -norm|n          Normalisation method (default=vsn)
  -species|s       Species name any standard e! species alias(will be reset to dbname/latin name e.g. "homo_sapiens")
  -location        Physical location of experimental group
  -contact         Contact details for experimental group
  -debug           Debug level (1-3)
  -log_file        Defines the log file
  -ucsc_coords     Flag to define usage of UCSC coord system in source files
  -release         Release version to load the registry from, defaults to latest release.
  -assembly        The specific assembly version to use e.g. 36 or 37.

  -update_xml      ?? Is this required, isn't this always updated now?


  -help            Brief help message
  -man             Full documentation


=head1 DESCRIPTION

B<This program> take several options, including an definitions file to parse and import array data into the ensembl-efg DB

=cut


#add @INC stuff here, or leave to .bashrc/.efg?

BEGIN{
	if(! defined $ENV{'EFG_DATA'}){
		if(-f "~/src/ensembl-functgenomics/scripts/.efg"){
			system (". ~/src/ensembl-functgenomics/scripts/.efg");
		}else{
			die ("This script requires the .efg file available from ensembl-functgenomics\n".
				 "Please source it before running this script\n");
		}
	}
}
	

#use Bio::EnsEMBL::Root; #Only used for rearrange see pdocs
#Roll own Root object to handle debug levels, logging, dumps etc.

### MODULES ###
use Getopt::Long;
use Pod::Usage;
use File::Path;
use Bio::EnsEMBL::Funcgen::Importer;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );

use strict;

$| = 1;#autoflush
my ($input_name, $input_dir, $name, $rset_name, $output_dir, $loc, $contact, $group, $pass, $dbname, $ssh);
my ($assm_ver, $help, $man, $species, $nmethod, $dnadb, $array_set, $array_name, $vendor, $exp_date, $ucsc);
my ($ctype, $ftype, $recover, $mage_tab, $update_xml, $write_mage, $no_mage, $farm, $exp_set, $old_dvd_format);
my ($reg_host, $reg_user, $reg_port, $reg_pass);
my ($parser, $fanal, $release, @result_files);

#to be removed
my ($import_dir);
my $reg = "Bio::EnsEMBL::Registry";
#?

my $data_dir = $ENV{'EFG_DATA'};
my $interactive = 1;
my $format = "tiled";
my $user = "ensadmin";
my $host = 'localhost';
my $port = '3306';
my $fasta = 0;#Shouldn't this be on by default?
my $verbose = 0;

#Definitely need some sort of Defs modules for each array?

$main::_debug_level = 0;
$main::_tee = 0;
my @tmp_args = @ARGV;
#Use some sort of DBDefs for now, but need  to integrate with Register, and have put SQL into (E)FGAdaptor?
#Use ArrayDefs.pm module for some of these, class, vendor, format?
#ArrayDefs would also contain paths to data and vendor specific parse methods?


GetOptions (
			#Experiment
			"name|n=s"           => \$name,
			"format|f=s"         => \$format,
			"vendor|v=s"         => \$vendor,
			"parser=s"           => \$parser,
			"array_set"          => \$array_set,
			"array_name=s"       => \$array_name,
			"result_set=s"       => \$rset_name,
			"experimental_set=s" => \$exp_set,
			"feature_type=s"     => \$ftype,
			"feature_analysis=s" => \$fanal,
			"cell_type=s"        => \$ctype,
			"exp_date=s"         => \$exp_date,
			'result_files=s{,}'  => \@result_files,

			#Experimental group 
			"group|g=s"    => \$group,
			"location=s"   => \$loc,
			"contact=s"    => \$contact,

			#Run mode
			"fasta"          => \$fasta,
			"farm=s"         => \$farm,
			"recover|r"      => \$recover,
			"norm_method=s"  => \$nmethod,
			"old_dvd_format" => \$old_dvd_format,
			"ucsc_coords"    => \$ucsc,
			"interactive"    => \$interactive,

			#MAGE
			"write_mage"   => \$write_mage,
			'update_xml'   => \$update_xml,
			'no_mage'      => \$no_mage,

			#DB connection
			"dbname=s"        => \$dbname,
			"pass|p=s"        => \$pass,
			"port|l=s"        => \$port,
			"host|h=s"        => \$host,
			"user|u=s"        => \$user,
			"species|s=s"     => \$species,
			"registry_user=s" => \$reg_user,
			"registry_host=s" => \$reg_host,
			"registry_pass=s" => \$reg_pass,
			"registry_port=s" => \$reg_port,
			"ssh"             => \$ssh,
			"assembly|a=s"    => \$assm_ver,
			"release=s"       => \$release,
		
			#Directory overrides
			"data_root=s"  => \$data_dir,
			"input_dir=s"  => \$input_dir,
			"output_dir=s" => \$output_dir,
		
			#Other params
			"tee"          => \$main::_tee,
			"log_file=s"   => \$main::_log_file,
			"debug_file=s" => \$main::_debug_file,
			"debug=i"      => \$main::_debug_level,		
			"help|?"       => \$help,
			"man|m"        => \$man,
			"verbose"      => \$verbose, # not implmented yet in Importer?
		   )   or pod2usage( -exitval => 1,
							 -message => "Params are:\t@tmp_args"
						   );

#Need to put these in an opt so we can test for the rest of ARGV i.e. we have missed an opt

#my @result_files = @ARGV;

if (@ARGV){
  pod2usage( -exitval =>1,
			 -message => "You have specified incomplete options:\t@tmp_args");
}

throw("Nimblegen import does not support cmdline defined result files") if (@result_files && uc($vendor) eq "NIMBELGEN");

#Need to add (primary) design_type and description, or add to defs file?


pod2usage(0) if $help;#Just synopsis
pod2usage(-exitstatus => 0, -verbose => 2) if $man;#All Pod with highlighting

#checking params used in this scripts
die('Must provide a -vendor parameter') if ! defined $vendor;
die('Must provide a -name parameter') if ! defined $name;
die('Must provide a -species parameter') if ! defined $species;


#log/debug files fail in Helper without this
$output_dir  = $data_dir."/output/${species}/".uc($vendor)."/".$name;
system("mkdir -p $output_dir -m 0755");
chmod 0755, $output_dir;


# pass as args?
$main::_log_file = $output_dir."/${name}.log" if(! defined $main::_log_file);
$main::_debug_file = $output_dir."/${name}.dbg" if(! defined $main::_debug_file);

### SET UP IMPORTER (FUNCGENDB/DNADB/EXPERIMENT) ###




my $Imp = Bio::EnsEMBL::Funcgen::Importer->new
  (
   -name        => $name,
   -format      => $format,
   -vendor      => $vendor,
   -parser      => $parser,
   -group       => $group,
   -pass        => $pass,
   -host        => $host,
   -user        => $user,
   -port        => $port,
   -registry_pass => $reg_pass,
   -registry_host => $reg_host,
   -registry_user => $reg_user,
   -registry_port => $reg_port,
   -ssh         =>  $ssh,
   -dbname      => $dbname,
   -array_set   => $array_set,
   -experimental_set_name => $exp_set,
   -array_name  => $array_name,
   -result_set_name => $rset_name, #not implemented yet
   -feature_type_name => $ftype,
   -feature_analysis => $fanal,
   -cell_type_name => $ctype,
   -write_mage    => $write_mage,
   -update_xml => $update_xml,
   -no_mage => $no_mage,
   -assembly => $assm_ver,
   -data_root   => $data_dir,
   -output_dir  => $output_dir,
   -recover     => $recover,
   -dump_fasta  => $fasta,
   -norm_method => $nmethod,
   -species     => $species,
   -farm        => $farm,
   -location    => $loc,
   -contact     => $contact,
   -verbose     => $verbose,
   -input_dir   => $input_dir,
   -exp_date     => $exp_date,
   -result_files => \@result_files,
   -old_dvd_format => $old_dvd_format,
   -ucsc_coords => $ucsc,
   -release => $release,
   #Exp does not build input dir, but could
   #This allows input dir to be somewhere 
   #other than efg dir structure
  );


$Imp->register_experiment();


1;
