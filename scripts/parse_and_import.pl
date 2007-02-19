#!/usr/local/ensembl/bin/perl -w

####!/opt/local/bin/perl -w


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
  -dbname          Defines the eFG dbname if it is not standard
  -ssh             Forces use of tcp protocol for using ssh forwarded ports from remote machine(remove and just specify 127.0.0.1 as host?)
  -array_set       Flag to treat all chip designs as part of same array
  -array_name      Name of the set of array chips
  -result_set      Name to give the raw/normalised result set.
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
#use Carp;#For dev only? cluck not exported by default Remove this and implement in Helper
use Pod::Usage;
#POSIX? File stuff
use File::Path;
use Bio::EnsEMBL::Funcgen::Importer;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );

use strict;

$| = 1;#autoflush
my ($input_name, $input_dir, $name, $rset_name, $output_dir, $loc, $contact, $group, $pass, $dbname, $ssh);
my ($data_version, $help, $man, $species, $nmethod, $dnadb, $array_set, $array_name, $vendor, $exp_date);
my ($ctype, $ftype, $recover);
my $reg = "Bio::EnsEMBL::Registry";

#to be removed
my ($import_dir);

my $data_dir = $ENV{'EFG_DATA'};
my $interactive = 1;
my $format = "tiled";
my $user = "ensadmin";
my $host = 'localhost';
my $port = '3306';
my $fasta = 0;#Shouldn't this be on by default?
#my $recover = 0;
my $verbose = 0;

#Definitely need some sort of Defs modules for each array?

$main::_debug_level = 0;
$main::_tee = 0;


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
	    "ssh"          => \$ssh,
	    "dbname=s"     => \$dbname,
	    "group|g=s"    => \$group,
	    "species|s=s"      => \$species, # Not needed as the db is species specific
	    "data_version|d=s" => \$data_version,
	    "array_set"    => \$array_set,
	    "array_name=s" => \$array_name,
	    "result_set=s" => \$rset_name,
	    "feature_type=s" => \$ftype,
	    "cell_type=s"    => \$ctype,
	    "debug=i"    => \$main::_debug_level,
	    "data_root=s"  => \$data_dir,
	    "input_dir=s"  => \$input_dir,
	    "exp_date=s"   => \$exp_date,
	    "fasta"        => \$fasta,
	    "recover|r"    => \$recover,
	    "norm=s"       => \$nmethod,
	    "location=s"   => \$loc,
	    "contact=s"    => \$contact,
	    "tee"          => \$main::_tee,
	    "log_file=s"   => \$main::_log_file,
	    "debug_file=s" => \$main::_debug_file,
	    #should have MAGE flag here? or would this be format?
	    "interactive"  => \$interactive,
	    "help|?"       => \$help,
	    "man|m"        => \$man,
	    "verbose"    => \$verbose, # not implmented yet
	   );


my @result_files = @ARGV;



throw("Nimblegen import does not support cmdline defined result files") if (@result_files && uc($vendor) eq "NIMBELGEN");

#Need to add (primary) design_type and description, or add to defs file?


pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;

#Build and validate all these in Experiment::new? We only need these for importing/analysing???
$output_dir  = $data_dir."/output/".uc($vendor)."/".$name;

#mkpath $output_dir;#log/debug files fail in Helper without this
system("mkdir -p $output_dir -m 0755");


chmod 0755, $output_dir;


# pass as args?
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
					       ssh         =>  $ssh,
					       dbname      => $dbname,
					       array_set   => $array_set,
					       array_name  => $array_name,
					       result_set_name => $rset_name, #not implemented yet
					       feature_type_name => $ftype,
					       cell_type_name => $ctype,
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
					       input_dir   => $input_dir,
					       exp_date     => $exp_date,
					       result_files => \@result_files,
					       #Exp does not build input dir, but could
					       #This allows input dir to be somewhere 
					       #other than efg dir structure
					      );



#Can be moved to set up script?

#print "exp is ".$Imp->db->get_ExperimentAdaptor->fetch_by_dbID(1)->name()."\n";

#exit;

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
				       -displayable     => 1,
				      );


my $sanger_anal = Bio::EnsEMBL::Analysis->new(
					      -logic_name      => 'SangerPCR',
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
					      -description    => 'Ratio generated by standard Sanger PCR array processing',
					      -display_label   => 'SangerPCR',
					      -displayable     => 1,
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
					     -displayable     => 1,
										  );

my $vsn_anal = Bio::EnsEMBL::Analysis->new(
										   -logic_name      => 'VSN_GLOG',
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
										   -description    => 'Generalised log transformation based on VSN variance stabilised scores',
										   -display_label   => 'VSN_GLOG',
					   -displayable     => 1,
										  );


my $t_anal = Bio::EnsEMBL::Analysis->new(
									   -logic_name      => 'Nessie',
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
									   -description    => 'Hidden Markov Model based predictions based on tiling array data',
									   -display_label   => 'Nessie (TilingHMM)',
					 -displayable     => 1,
									  );
my $lanal = Bio::EnsEMBL::Analysis->new(
				       -logic_name      => 'LiftOver',
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
				       -description    => 'Remapping to new assembly performed by LiftOver',
				       -display_label   => 'LiftOver',
				       -displayable     => 1,
				      );


#this checks if already stored
$anal_a->store($anal);
$anal_a->store($lanal);
$anal_a->store($raw_anal);
$anal_a->store($vsn_anal);
$anal_a->store($t_anal);
$anal_a->store($sanger_anal);


#Validate, parse and import all experiment data


$Imp->register_experiment();


1;
