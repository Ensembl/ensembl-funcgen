#!/usr/bin/env perl

=head1 LICENSE


  Copyright (c) 1999-2011 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

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
                     This will be used to select the import parsers if --parser is not defined.
  --parser           Selects the import parser to use e.g. Bed etc (Will override default vendor --parser)
  --array_name       Name of the set of array chips(Also see -array_set option)
  --result_set       Name to give the raw/normalised result set.
  --feature_type     The name of the FeatureType of the experiment e.g. H4K36me3.
  --input_set        Defines the name of an InputSet import
  --input_feature_class Defines the type of features being loaded via an InputSet(to be used 
                     with -input_set e.g. result or annotated).
  --_total_features  Specifies the total numbers of feature in the InputSet for calculating RPKM.
                     This is normally handled by the prepare stage and never has to specified by the user.
  --cell_type        The name of the CellType of the experiment.
  --feature_analysis The name of the Analysis used in the experiment.
  --norm|n           Normalisation method (default=vsn)
                     NOTE: FeatureType, CellType and Analysis(inc norm) entries must already 
                     exist in the eFG DB. See ensembl-functgenomics/script/import/import_types.pl
  --exp_date         The date for the experiment.
  --result_files     Space separated list of result files paths (Used in InputSet imports e.g. Bed).
  --slices           Space separated list of slice names|seq_region_names to import (only works for -input_set import)
  --skip_slices      Space separated list of seq_region names to skip (only works for -input_set import)

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
  --queue           LSF queue to submit farm jobs to
  --no_log                 
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
  --lsf_name        Name of host as monitored by lsf, default is host =~ s/-/_/ and prefixed with 'my' e.g. myensgenomics_1

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
  --no_log          Prevents log file form being written, used to capture log in lsf out file.
  --help            Brief help message
  --man             Full documentation
  --verbose


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
	

### MODULES ###
use Getopt::Long;
use Pod::Usage;
use File::Path;
use Bio::EnsEMBL::Funcgen::Importer;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw (strip_param_args generate_slices_from_names strip_param_flags);
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw( warning );
use strict;

$| = 1;#autoflush
my ($input_name, $input_dir, $name, $rset_name, $output_dir, $loc, $contact, $group, $pass, $dbname, $ssh);
my ($assm_ver, $help, $man, $species, $nmethod, $dnadb, $array_set, $array_name, $vendor, $exp_date, $ucsc);
my ($ctype, $ftype, $recover, $mage_tab, $update_xml, $write_mage, $no_mage, $farm, $exp_set, $old_dvd_format);
my ($reg_host, $reg_user, $reg_port, $reg_pass, $input_feature_class, $lsf_host, $batch_job, $prepared, $total_features);
my ($parser, $fanal, $release, $format, $on_farm, @result_files, @slices, @skip_slices);

my $data_dir = $ENV{'EFG_DATA'};
my $interactive = 1;
my $user = $ENV{'DB_USER'};
my $host = 'localhost';
my $port = '3306';
my $fasta = 0;
my $verbose = 0;
my $queue = 'long';

#Helper config
$main::_debug_level = 0;
$main::_tee = 0;


my @tmp_args = @ARGV;


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
			'skip_slices=s{,}'   => \@skip_slices,
			'_total_features=i'   => \$total_features,

			#Experimental group 
			"group|g=s"    => \$group,
			"location=s"   => \$loc,
			"contact=s"    => \$contact,

			#Run modes
			'fasta'          => \$fasta,
			'farm'           => \$farm,
			#Internal run modes
			'_batch_job'      => \$batch_job,
			'_on_farm'        => \$on_farm,
			'_prepared'       => \$prepared, 
			#Internal option to signify job has been previously prepared
			#i.e. result_file name may differ from the original InputSubset name
			#slightly different to batch_job as prepare may not always result in a changed file
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
			"lsf_host=s"      => \$lsf_host,
			"registry_user=s" => \$reg_user,
			"registry_host=s" => \$reg_host,
			"registry_pass=s" => \$reg_pass,
			"registry_port=s" => \$reg_port,
			"ssh"             => \$ssh,
			"assembly|a=s"    => \$assm_ver,
			"release=s"       => \$release,
		
			#Directory overrides
			"data_root=s"       => \$data_dir,
			"input_dir=s"       => \$input_dir,
			"output_dir=s"      => \$output_dir,

			#Other params
			"help|?"       => \$help,
			"man|m"        => \$man,
			"verbose"      => \$verbose, # not implmented yet in Importer?
			#Helper config
			"tee"          => \$main::_tee,
			"log_file=s"   => \$main::_log_file,
			"no_log"       => \$main::_no_log,
			"debug_file=s" => \$main::_debug_file,
			"debug=i"      => \$main::_debug_level,		
		   )
  or pod2usage( -exitval => 1,
				-message => "Params are:\t@tmp_args"
			  );


print "parse_and_import.pl @tmp_args\n";


if($batch_job && ! defined $ENV{'LSB_JOBINDEX'}){
  die('Your -_batch_job is not running on the farm. Did you try and set this internal parameter manually?');
}


if (@ARGV){
  pod2usage( -exitval =>1,
			 -message => "You have specified unknown options/args:\t@ARGV");
}

die("Nimblegen import does not support cmdline defined result files") if (@result_files && uc($vendor) eq "NIMBELGEN");

#Need to add (primary) design_type and description, or add to defs file?


pod2usage(0) if $help;#Just synopsis
pod2usage(-exitstatus => 0, -verbose => 2) if $man;#All Pod with highlighting

#checking params used in this scripts
die('Must provide a -vendor parameter') if ! defined $vendor;
die('Must provide a -name parameter') if ! defined $name;
die('Must provide a -species parameter') if ! defined $species;


#log/debug files fail in Helper without this
$output_dir  ||= $data_dir."/output/${dbname}/".uc($vendor)."/".$name;
system("mkdir -p $output_dir -m 0755");
chmod 0755, $output_dir;

if(! $main::_no_log){
  $main::_log_file = $output_dir."/${name}.log" if(! defined $main::_log_file);
  $main::_debug_file = $output_dir."/${name}.dbg" if(! defined $main::_debug_file);
}


#Define import parser
#This will override the default Vendor Parser type
#Evals simply protect from messy errors if parser type not found
my $parser_error;
my $vendor_parser =  ucfirst(lc($vendor));

eval {require "Bio/EnsEMBL/Funcgen/Parsers/${vendor_parser}.pm";};
 
if($@){
  #Don't warn/throw yet as we might have a standard parser format
  
  $parser_error .= "There is no valid parser for the vendor your have specified:\t".$vendor.
	"\nMaybe this is a typo or maybe you want to specify a default import format using the -parser option\n".$@;
}

if(defined $parser){
  
  #try normal case first
  eval {require "Bio/EnsEMBL/Funcgen/Parsers/${parser}.pm";};
  
  if($@){
	$parser = ucfirst(lc($parser));
	
	#Now eval the new parser
	eval {require "Bio/EnsEMBL/Funcgen/Parsers/${parser}.pm";};
	
	if($@){
	  my $txt = "There is no valid parser for the -parser format your have specified:\t".$parser."\n";
	  
	  if(! $parser_error){
		$txt .= "Maybe this is a typo or maybe you want run with the default $vendor_parser parser\n";
	  }
		
	  die($txt.$@);
	}
	
	#warn about over riding vendor parser here
	if(! $parser_error){
	  #Can't log this as we haven't blessed the Helper yet
	  warn("WARNING\t::\tYou are over-riding the default ".$vendor." parser with -parser ".$parser);
	}
  }
}
else{
  die($parser_error) if $parser_error;
  $parser = $vendor_parser;
}


if($input_feature_class eq 'result' && 
   ((! @slices) || (scalar(@slices) > 1)) && #have more than one slice
   (! ($farm || $batch_job))){
  die('ResultFeature Collections cannot be imported across all slices in one job.  It is more sensible to submit parallel jobs to the farm using the -farm option'); 
}

if(@slices && ! $exp_set){
  die('-slices parameter is only valid for -input_set import');
}


### SET UP IMPORTER (FUNCGENDB/DNADB/EXPERIMENT) ###

#Only change this to $parser->new when we have implemented the changes in the parsers.
#This involves changing the inheritance
#and moving the set_defs methods into the new method

my $Imp = Bio::EnsEMBL::Funcgen::Importer->new
  (
   -name        => $name,
   -format      => $format,
   -vendor      => $vendor,
   -parser      => $parser,##########################REMOVE THIS WHEN WE HAVE MADE THE CHANGES IN THE PARSER INHERITANCE
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
   -batch_job   => $batch_job,
   -prepared    => $prepared,
   -location    => $loc,
   -contact     => $contact,
   -verbose     => $verbose,
   -input_dir   => $input_dir,
   -exp_date     => $exp_date,
   -result_files    => \@result_files,
   -total_features => $total_features,
   -old_dvd_format => $old_dvd_format,
   -ucsc_coords => $ucsc,
   -release => $release,
   #Exp does not build input dir, but could
   #This allows input dir to be somewhere 
   #other than efg dir structure
  );

print "The log files and output can be found here:\t".$Imp->get_dir('output')."\n";


#Need to think about how best to handle this job submission
#Convert to pipeline?
#-top_level_slices
#-farm
#-no_farm

#do not force farm here as people may not have access

my $slice_adaptor = $Imp->slice_adaptor;

if(@slices || $input_feature_class eq 'result'){

  if($input_feature_class ne 'result'){
	warn "You are setting slices for an input_feature_class which does not support slice based import:\t$input_feature_class\n";
  }

  if(! @slices){
	$Imp->log("No slices defined defaulting to current toplevel");
  }

  @slices = @{&generate_slices_from_names($slice_adaptor, \@slices, \@skip_slices, 1, undef, 1)};#toplevel, nonref, incdups
  #inc dups here for now for loading Y
  #Until we support PAR/HAP regions properly
}



### PREPARE

my $input_file  = '';
my $total_feats = '';

# Only need to prepare for result input_sets
# which have not been previously prepared
# We always need to prepare now as RPKM requires total_feature count

if((defined $input_feature_class) &&
   ($input_feature_class eq 'result') &&
   (! $prepared)){
  
  #Sanity checks
  if($batch_job){
	die('You should prepare prior to submitting a -_batch_job. Have you set this private param manually?');
  }

  if($total_features){
	die("-total_feature is already specified($total_features), did you set this manually? This should be set by the 'prepare' step.");
  }

    
  $Imp->slices(\@slices);       #Set slices so we can preprocess the file
  $Imp->init_experiment_import; #Define and store sets once here before setting off parallel farm jobs.

  print "Preparing data...\n";
  $Imp->read_and_import_data('prepare');
  $prepared = ' -_prepared ';
  #This prevents 'No space left on device' errors about /tmp space
  #Parses, filters and sorts input to tmp file
  #Slight overhead in certain circumstances, But is cleaner and more 
  #likely to fail here instead of in mutiple batch jobs
  #Could output to seq_region named files which would speed up 2nd stage parsing
  #Delete immediately after use, or save and md5?
  #Do we need to put these in a large file stripe dir?
  #Same with mapping pipeline?


  #Now reset the slices based on what we have seen in the data
  #This will prevent 'empty' collection blobs
  @slices = values(%{$Imp->slice_cache});
  $Imp->slices(\@slices);

  if(! ($input_file = $Imp->output_file)){

	#Die for now as we currently always expect a sorted file
	#and we are not catching some failures here
	#Is this to do with?
	#gzip: stdout: Disk quota exceeded

	die("Could not get prepared file for:\t".$result_files[0]);

	#$input_file = $result_files[0];
	#$prepared = '';
  }

  #Grab total features count as prepared file may have been filtered for specified slices
  $total_feats = $Imp->counts('total_features');

  if(! $total_feats){
	die("Could not get 'total_features' from InputSet Parser");
  }

  #Re/set attrs/params appropriately for local/farm job
  if($farm){
	$total_feats = "-_total_features $total_feats";
	$input_file = "-result_files $input_file";
  }
  else{
	#re/set Imp attrs if we are running locally
	$Imp->total_features($total_feats);
	$Imp->prepared(1);
	$Imp->result_files([$input_file]);
  }
}
else{
  #Need to format vars as above

  #Logically these may never be needed in the farm submission
  #but define in case we change the logic
  $total_feats = "-_total_features $total_features" if $total_features;
  $input_file  = "-result_files @result_files" if @result_files;
  
}


#Can we change the logic of these tests to allow passing of _total_features and _prepared
#to avoid prepare step?
#This would also require regenerating the file name

### SUBMIT/RUN JOB(S) 

if($farm &&           ###BSUB JOBS
   (! $batch_job) &&
   (! $on_farm)){
  
  #We don't strip farm from args, so this will currently loop
  #we use farm and batch_job in importer
  my $slices    = '';
  my $batch_job = '';
  
  if(@slices){
	$slices    = '-slices '.join(' ', (map $_->seq_region_name, @slices));
	$batch_job = '-_batch_job';
  }


  #STRIP ARGS
  #strip slice params as we have already defined which slices to use
  #strip log/debug_flie and set -no_log as we will use the lsf output
  #redefine result file as prepared file
  #set new slices as those which have been seen in the input
  #so we don't have to run jobs for all toplevel seq_regions
  my @args = @{&strip_param_args(\@tmp_args, ('log_file', 'debug_file', 
											  'result_files', 'slices', 
											  'skip_slices'))};
  #no log to avoid overwriting log file
  #output in lsf files
  #should probably change this by setting log dependant on batch_job in Importer

  my $cmd = "time perl $ENV{EFG_SRC}/scripts/import/parse_and_import.pl -_on_farm -no_log $batch_job @args $input_file $prepared $total_feats $slices";

  #sub slices for keys in slice cache? This will lose any incomplete slice names?
  #But we don't allow partial Slice imports yet!
  #So it is fine to use the slice_cache keys
  



  #Set up batch job info
  #Batches obnubilate the job info/output slightly as we don't know which slice is running

  my $job_name = "parse_and_import_${name}";

  #THROTTLE
  #Load will only ever be 12 max, 1 for each connection and 10 for an active query
  #Dropped duration to 5 as this only sets the interval at which LSF tries to submit jobs
  #before using the select load level to submit any further jobs
  #Is the load select too high or duration too low? We are getting a total load of 800-900 with ~130 jobs running
  #Dropped to 600.

  if(! $lsf_host){
	($lsf_host = 'my'.$host) =~ s/-/_/g;
  }

  

  my $rusage = "-R \"select[($lsf_host<=600)] rusage[${lsf_host}=12:duration=10]\" ".
	"-R \"select[mem>6000] rusage[mem=6000]\" -M 6000000";
  #Upped the duration as the first stage is to preprocess the file and will not show load on the server


  my $bsub_cmd;

  if(@slices){
	#limit concurrent batch jobs, added as lsf throttling not working
   	#warn "HARDCODING job limit a 2/array";
	#$bsub_cmd="bsub $rusage -q $queue -J \"${job_name}[1-".scalar(@slices)."]%2\" -o ${output_dir}/${job_name}.".'%J.%I'.".out -e ${output_dir}/${job_name}.".'%J.%I'.".err $cmd";
	$bsub_cmd="bsub $rusage -q $queue -J \"${job_name}[1-".scalar(@slices)."]\" -o ${output_dir}/${job_name}.".'%J.%I'.".out -e ${output_dir}/${job_name}.".'%J.%I'.".err $cmd";
	
  }
  else{
	$bsub_cmd="bsub $rusage -q $queue -J \"${job_name}\" -o ${output_dir}/${job_name}.out -e ${output_dir}/${job_name}.err $cmd";
  }


  $Imp->log("Submitting $job_name\n$bsub_cmd");
  $Imp->log("Slices:\n\t".join("\n\t", (map $_->name, @slices))."\n") if @slices;
  system($bsub_cmd) && die("Failed to submit job:\t$job_name\n$bsub_cmd");

   
  warn "Can't set imported states for InputSubset and ResultSet until jobs have completed succesfully, please update manually for now.";
  #Imported state may also be present even when we set a subset of slices!
  #Hence wwe won't be able to import more slice data from the same file if we set it as IMPORTED!
  #We currently get around this by never setting it, and always rolling back Slice based imports

  #Wait here before:
  #   checking job success
  #   setting any extra slice specific status entries


}
else{ ### DO THE IMPORT
  print "Importing data...\n";

  @slices = ($slices[($ENV{LSB_JOBINDEX} - 1)]) if $batch_job;

  #warn "LSB_JOBINDEX is ".$ENV{LSB_JOBINDEX};

  #warn "New slices are:\t@slices";

  $Imp->slices(\@slices) if @slices;

  #warn "Importer Slices are ".join(', ', @{$Imp->slices(\@slices)});

  $Imp->register_experiment();
}

print "Done\n";

1;
