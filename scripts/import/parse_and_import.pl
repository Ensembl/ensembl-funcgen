#!/software/bin/perl -w


=head1 NAME

 ensembl-efg parse_and_import.pl
  
=head1 SYNOPSIS

 parse_and_import.pl [options]

=head1 OPTIONS

 Options:

 Experiment (Mostly Mandatory)
  --name|n           Instance/Experiment name. This will be used as the input/out_dir 
                     unless otherwise specified (--input_dir/--output_dir).
  --format|f         Format of experiment technology e.g. TILING.
  --vendor|v         Vendor of experiment technology e.g. NIMBLEGEN
                     This will be used to select the import parsers if the --parser
                     is not defined.
  --parser           Selects the import parser to use e.g. NIMBLEGEN, Bed, Sanger.
  --array_name       Name of the set of array chips(Also see -array_set option)
  --result_set       Name to give the raw/normalised result set.
  --experimental_set Name to give the ExperimentalSet for 
  --feature_type     The name of the FeatureType of the experiment e.g. H4K36me3.
  --cell_type        The name of the CellType of the experiment.
  --feature_analysis The name of the Analysis used in the experiment.
  --norm|n           Normalisation method (default=vsn)
                     NOTE: FeatureType, CellType and Analysis(inc norm) entries must already 
                     exist in the eFG DB. See ensembl-functgenomics/script/import/import_types.pl
  --exp_date         The date for the experiment.
  --result_files     Space separated list of result files paths (Used in Bed import).

 Experimental group (Mostly Mandatory)
  --format|f        Data format
  --group|g         Group name
  --location        Physical location of experimental group
  --contact         Contact details for experimental group

 Run modes
  --array_set       Flag to treat all chip designs as part of same array
  --ucsc_coords     Flag to define usage of UCSC coord system in source files (chr_start==0).
  --fasta           Fasta dump flag
  --farm
  --interactive
  --old_dvd_format  Flag to use the old NIMBLEGEN DVD format
  --recover         Flag to enable rollback and over-writing of previously imported data
                    NOTE: This must be on to run the 2nd stage of a MAGE inport e.g. NIMBLEGEN


 MAGE
  ##--update_xml      Deprecated, no longer used as we always update now
  --write_mage      Flag to force 1st stage of MAGE import i.e. writing of 
                    tab2mage file to define replicates and experimental meta data.
  --no_mage         Flag to run import without mage creation/validation
                    WARNING: Not advised as this can cause erroneous replicate ResultSet generation

 eFG DB Connection (Mostly Mandatory)
  --dbname          The eFG dbname
  --user            The eFG user
  --pass|p          The eFG write password 
  --port            The eFG port
  --host            The eFG host
  --species|s       Species name any standard e! species alias(will be reset to dbname/latin name e.g. "homo_sapiens")
 
  --ssh             Forces use of tcp protocol for using ssh forwarded ports from remote machine(remove and just specify 127.0.0.1 as host?)

 Core DB Connection (This will default to load from ensembldb.ensembl.org)
  --registry_host   Host to load registry from
  --registry_port   Port for registry host
  --registry_user   User for registry host 
  --registry_pass   Password for registry user
  --release         Release version to load the registry from, defaults to latest release.
  --assembly        The specific assembly version to use e.g. 36 or 37.
 

 Directory over-rides
  --data_root       The root data dir ($ENV{'EFG_DATA'}).
  --input_dir       Over-rides use of -name to define the input path
  --output_dir      Over-rides use of -name to define the output path
  
 Other
  --tee             Outputs logs to STDOUT aswell as log file.
  #--debug           Debug level (1-3). Not implemented
  --log_file        Defines the log file, default is $output_dir/"epxeriment_name".log
  #--debug_file      Defines the debug file. Not implemented
  --help            Brief help message
  --man             Full documentation
  --verbose


=head1 DESCRIPTION

B<This program> take several options, including an definitions file to parse and import array data into the ensembl-efg DB


=head1 LICENSE

  Copyright (c) 1999-2009 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.


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
			"array_name=s"       => \$array_name,
			"result_set=s"       => \$rset_name,
			"experimental_set=s" => \$exp_set,
			"feature_type=s"     => \$ftype,
			"feature_analysis=s" => \$fanal,
			"cell_type=s"        => \$ctype,
			"norm_method=s"      => \$nmethod,
			"exp_date=s"         => \$exp_date,
			'result_files=s{,}'  => \@result_files,

			#Experimental group 
			"group|g=s"    => \$group,
			"location=s"   => \$loc,
			"contact=s"    => \$contact,

			#Run modes
			"fasta"          => \$fasta,
			"farm=s"         => \$farm,
			"recover|r"      => \$recover,
			"old_dvd_format" => \$old_dvd_format,
			"ucsc_coords"    => \$ucsc,
			"interactive"    => \$interactive,
			"array_set"          => \$array_set,

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
