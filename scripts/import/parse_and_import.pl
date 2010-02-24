#!/software/bin/perl -w


=head1 NAME

 ensembl-efg parse_and_import.pl
  
=head1 SYNOPSIS

 parse_and_import.pl [options]

=head1 OPTIONS

 Options:

 Experiment/Set (Mostly Mandatory)
  --name|n           Instance/Experiment name. This will be used as the input/out_dir if applicable,
                     unless otherwise specified (--input_dir/--output_dir).
  --format|f         Format of experiment technology e.g. TILING.
  --vendor|v         Vendor of experiment technology e.g. NIMBLEGEN
                     This will be used to select the import parsers if the --parser
                     is not defined.
  --parser           Selects the import parser to use e.g. NIMBLEGEN, Bed, Sanger.
  --array_name       Name of the set of array chips(Also see -array_set option)
  --result_set       Name to give the raw/normalised result set.
  --feature_type     The name of the FeatureType of the experiment e.g. H4K36me3.
  --input_set        Defines the name of an InputSet import
  --input_feature_class Defines the type of features being loaded via an InputSet(to be used 
                     with -input_set e.g. result or annotated).
  --cell_type        The name of the CellType of the experiment.
  --feature_analysis The name of the Analysis used in the experiment.
  --norm|n           Normalisation method (default=vsn)
                     NOTE: FeatureType, CellType and Analysis(inc norm) entries must already 
                     exist in the eFG DB. See ensembl-functgenomics/script/import/import_types.pl
  --exp_date         The date for the experiment.
  --result_files     Space separated list of result files paths (Used in InputSet imports e.g. Bed).
  --slices           Space separated list of slice names to import(Only works for -input_set import)

 Experimental group (Mostly Mandatory)
  --format|f        Data format
  --group|g         Group name
  --location        Physical location of experimental group
  --contact         Contact details for experimental group

 Run modes
  --array_set       Flag to treat all chip designs as part of same array
  --ucsc_coords     Flag to define usage of UCSC coord system in source files (chr_start==0).
  --fasta           Fasta dump flag
  --farm            Dependant on the import type, this either submits slice based import jobs or 
                    array normalisation jobs to the farm
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
  --log_file        Defines the log file, default is $output_dir/"epxeriment_name".log
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
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw (strip_param_args generate_slices_from_names);
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );

use strict;

$| = 1;#autoflush
my ($input_name, $input_dir, $name, $rset_name, $output_dir, $loc, $contact, $group, $pass, $dbname, $ssh);
my ($assm_ver, $help, $man, $species, $nmethod, $dnadb, $array_set, $array_name, $vendor, $exp_date, $ucsc);
my ($ctype, $ftype, $recover, $mage_tab, $update_xml, $write_mage, $no_mage, $farm, $exp_set, $old_dvd_format);
my ($reg_host, $reg_user, $reg_port, $reg_pass, $input_feature_class);
my ($parser, $fanal, $release, @result_files, @slices);

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
my $queue = 'long';

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
			"input_set=s"        => \$exp_set,
			"input_feature_class=s" => \$input_feature_class,
			"feature_type=s"     => \$ftype,
			"feature_analysis=s" => \$fanal,
			"cell_type=s"        => \$ctype,
			"norm_method=s"      => \$nmethod,
			"exp_date=s"         => \$exp_date,
			'result_files=s{,}'  => \@result_files,
			'slices=s{,}'        => \@slices,


			#Experimental group 
			"group|g=s"    => \$group,
			"location=s"   => \$loc,
			"contact=s"    => \$contact,

			#Run modes
			"fasta"          => \$fasta,
			"farm"           => \$farm,
			"queue=s"        => \$queue,
			"recover|r"      => \$recover,
			"old_dvd_format" => \$old_dvd_format,
			"ucsc_coords"    => \$ucsc,
			"interactive"    => \$interactive,
			"array_set"      => \$array_set,

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
			 -message => "You have specified unknown options/args:\t@ARGV");
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

#If we have slices we need to submit these to the farm as individual jobs
#We need to define the DBs here before we can generate the slices
#Or can we just do the job submission in the InputSet?
#No we really need all the args here.

#Let's just deal with one slice for now.

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
   -input_set_name => $exp_set,
   -input_feature_class => $input_feature_class,
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


#Need to think about how best to handle this job submission
#Convert to pipeline?
#Also need skip slices here
#-top_level_slices
#-farm
#-no_farm

#do not force farm here as people may not have access

if($input_feature_class eq 'result' && 
   ((! @slices) || (scalar(@slices) > 1)) 
   && ! $farm){
  die('ResultFeature Collections cannot be imported across all slices in one job.  It is more sensible to submit parallel jobs to the farm using the -farm option'); 
}

my $slice_adaptor = $Imp->slice_adaptor;



if(@slices || $input_feature_class eq 'result'){

  if($input_feature_class ne 'result'){
	warn "You are setting slices for an input_feature_class which does not support slice based import:\t$input_feature_class\n";
  }

  if(! @slices){
	print "No slices defined defaulting to current toplevel\n";
  }

  @slices = @{&generate_slices_from_names($slice_adaptor, \@slices, 1, 1)};
}

#farm behaves differently if slices not defined
#i.e. submit norm jobs to farm from within Importer
#Will have to change this if we ever want to submit slice based 
#jobs from within the Importer

if(@slices && $farm){
  #submit slice jobs to farm
  #use Importer params to submit rather than @ARGV?
  #Strip -farm and restrict to just one slice

  #Can we put this in EFGUtils::bsub_by_Slice?
  #Or just remove args?

  my @args = @{&strip_param_args(\@tmp_args, ('slices'))};
  my $args = "@args";
  $args =~ s/ -farm( |$)//;
  my $cmd = "perl $ENV{EFG_SRC}/scripts/import/parse_and_import.pl $args";
  #qmy $lsf_out = $out_dir

  foreach my $slice(@slices){
	my $job_name = "parse_and_import_${name}_".$slice->name;

	#-R "select[mem>20000] rusage[mem=20000] -M 20000000"? 
	#Should probably use job arrays here as we already know how many jobs we are submitting, so we could use
	#$LSB_JOBINDEX to pick out the correct slice from an ordered array

	#Also need to redefine the log file as we will have many jobs writing to the same log file?
	#Or can we just turn logging off and write to lsf outfile?
	my $bsub_cmd="bsub -q $queue -J $job_name -o ${output_dir}/${job_name}.out -e ${output_dir}/${job_name}.err $cmd -slices ".$slice->name;

	print "Submitting $job_name\n";
	system($bsub_cmd) && die("Failed to submit job:\t$job_name\n$bsub_cmd");
  }
}
else{
  $Imp->slices(\@slices) if @slices;
  $Imp->register_experiment();
}


1;
