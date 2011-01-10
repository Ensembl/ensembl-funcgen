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

configure_inputs.pl - Generates a list of input_ids, configuring links
                      to input files where appropriate.


=head1 SYNOPSIS

configure_inputs.pl [ OPTIONS ]

=head1 OPTIONS

    -dbhost      
    -dbport      
    -dbuser      
    -dbpass      
    -dbname      
    -slice       Sets slice as input id type (encode|toplevel)
    -array       Sets array as input id type
    -dir         Sets file as input id type using files in dir input
    -logic_names List of analysis logic names e.g. SWEmbl ...
    -exp_regex   regular expression to select certain experiments
                      (default: fetch all available if omitted)
    -exp_suffix  
    -zip         Flag to zip input files if required

=head1 DESCRIPTION

This script generates a list of input_ids from Encode regions / files for
setting up the eFG peak calling pipeline.

=cut

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

my ($pass, $dnadb_pass, $species, $data_version,
	$exp_suffix, $slice_type, $array, $dir, $help, $man, $debug, $workdir,
	$host, $port, $user, $dbname, $dnadb_host, $zip,
	$dnadb_port, $dnadb_user, $dnadb_name, @logic_names);

my $exp_regex = '.*';


#warn "@ARGV";

#To do
# 1 Remove requirement fot gzip'd files. Just warn and zip them in situ or in the cache if they aren't already?
#   Need to add zip step to end of run_bwa.sh/read mapping pipeline.
# 2 Throw is not safe in script? Use die instead
# 3 Add $PDB_SCRIPT_ARGS(update these in pipeline env to -pdb_name etc?)
# 4 Better Pod and inline docs
# 5a DONE Move is_ methods and encode methods to EFGUtils, implement in SWEmbl Runnable for setting input file type
# 5b Polish above EFGUtils methods and write is_bam(needs adding to SWEmbl too)
# 6 Deprecate or update non-file based methods as these are currently untested.
# 7 Handle experiment/feature/data set naming better. Experiment name should be user defined or
#   parsed from dir in env. Set names should be parsed from CELL_TYPE_FEAUTURE_TYPE file name.
#   This is what exp_suffix is for, but we need to tidy this up i.e set exp_suffix by default in env.
# 8 Change suffix to prefix?

GetOptions (
            'dbpass|p=s'       => \$pass,
            'dbport=i'         => \$port,
            'dbhost|h=s'       => \$host,
            'dbuser|u=s'       => \$user,
            'dbname|d=s'       => \$dbname,
			'dnadb_pass=s'     => \$dnadb_pass,
            'dnadb_port=i'     => \$dnadb_port,
            'dnadb_host=s'     => \$dnadb_host,
            'dnadb_user=s'     => \$dnadb_user,
            'dnadb_name=s'     => \$dnadb_name,
            'species=s'        => \$species,
			'logic_names=s{,}' => \@logic_names,
            'data_version=s'   => \$data_version,
            'exp_regex|e=s'    => \$exp_regex,
            'exp_suffix=s'     => \$exp_suffix,
            'help|?'           => \$help,
            'man|m'            => \$man,
            'debug'            => \$debug,
            'slice=s'          => \$slice_type,
            'array'            => \$array,
            'dir=s'           => \$dir,
			'work_dir=s'       => \$workdir,
			'zip'              => \$zip,
			

           );

### defaults ###
if (!$port) {
	$port = 3306;
	warn("No port specified, using default '$port'.")
}
if (!$species) {
	$species = 'homo_sapiens';
	warn("No species specified, using default '$species'.")
}

### check options ###
die('Must provide a -logic_name') if ! @logic_names;
die('Must provide a -work_name') if (! $workdir || ! -d $workdir);
throw("Must specify mandatory database hostname (-host).\n") if ! defined $host;
throw("Must specify mandatory database username. (-user)\n") if ! defined $user;
throw("Must specify mandatory database name (-dbname).\n") if ! defined $dbname;
throw("Must specify mandatory password (-pass).\n") if ! defined $pass;

$| = 1;
#Need to remove use of throw here?
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info stack_trace_dump);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(open_file get_file_format is_gzipped);
use Bio::EnsEMBL::Funcgen::Utils::Encode qw(get_encode_regions);#move to EFGUtils

# Get eFG database adaptor

#This scripts should work independantly of the peak env
#remove these env vars
my $dnadb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (
     -host    => $dnadb_host,
     -user    => $dnadb_user,
     -port    => $dnadb_port,
     -dbname  => $dnadb_name,
     -species => $species,
     -pass    => $dnadb_pass,
     );


my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
    (
     -host    => $host,
     -user    => $user,
     -dbname  => $dbname,
     -species => $species,
     -pass    => $pass,
     -port    => $port,
     -dnadb   => $dnadb,
     );

#Need to remove this dependancy on env var
#print to file or add opts


my $pdb = Bio::EnsEMBL::Pipeline::DBSQL::DBAdaptor->new
    (
     -host    => $host,
     -user    => $user,
     -dbname  => $ENV{'PDB_NAME'},
     -species => $species,
     -pass    => $pass,
     -port    => $port,
     );




if ($slice_type) {
        
    # Get all experiments
    my $ea = $db->get_ExperimentAdaptor();
    my $exp = $ea->fetch_all();
        
    my @input_ids = ();

    if ($slice_type eq 'encode') {
        
        # Get Encode regions
        my $encode_regions = get_encode_regions($db->dnadb());
        #map { print Dumper $_->slice->name } @$encode_regions;
        
        foreach my $e (@$exp) {
            
            next if (defined $exp_regex && $e->name !~ m/$exp_regex/);
            map {
                push ( @input_ids, join(':', $e->name, $_->slice->coord_system->name, 
                                        $_->slice->coord_system->version,
                                        $_->slice->seq_region_name, $_->start,
                                        $_->end, $_->strand) );
            } @$encode_regions;
                
        }

    } elsif ($slice_type eq 'toplevel') {

        # Get toplevel slices
        my $sa = $db->get_SliceAdaptor();
        my $tls = $sa->fetch_all('toplevel');

        foreach my $e (@$exp) {

            next if (defined $exp_regex && $e->name !~ m/$exp_regex/);
            map {
                push ( @input_ids, join(':', $e->name, $_->name ) );
            } @$tls;
        }

    } else {
        throw("$slice_type is not a valid slice type. You need to specify either 'toplevel' or 'encode'.");
        
    }

    foreach my $input_id (@input_ids) {

	  foreach my $lname(@logic_names){
		
		my $analysis_id = &get_analysis_id($lname);

        my $sql = "insert into input_id_analysis (input_id,analysis_id,input_id_type)".
		  " values ('${input_id}',${analysis_id},'slice');";
        #warn($sql);
        eval {
            $pdb->dbc->do($sql);
		  };
        throw("Couldn't store input_id '$input_id'. Most likely it has already been ".
              "stored. Drop your input_ids with CleanInputIds and rerun CreateInputIds.")
		  if ($@);
	  }
	}


} elsif ($array) {

    # Get all experiments
    my $ea = $db->get_ExperimentAdaptor();
    my $exp = $ea->fetch_all();

    my @input_ids = ();

    foreach my $e (@$exp) {
        
        next if (defined $exp_regex && $e->name !~ m/$exp_regex/);
        
        #warn $e->name;
        push ( @input_ids, join(':', $e->name, "ARRAY" ) );

    }

    foreach my $input_id (@input_ids) {

	  foreach my $lname(@logic_names){
		my $analysis_id = &get_analysis_id($lname);

        my $sql = "insert into input_id_analysis (input_id,analysis_id,input_id_type)".
            " values ('${input_id}',${analysis_id},'array');";
        #warn($sql);
        eval {
		  $pdb->dbc->do($sql);
        };
        throw("Couldn't store input_id '$input_id'. Most likely it has already been ".
              "stored. Drop your input_ids with CleanInputIds and rerun CreateInputIds.")
		  if ($@);
	  }
	}

} elsif ($dir) {

	# get cell and feature type adapter
	my $cta = $db->get_CellTypeAdaptor();
	my $fta = $db->get_FeatureTypeAdaptor();

       
    if (! $dir) {
        throw("Need to specify a input directory containing ". 
              "files to be processed")
    }

    opendir(DIR, $dir)
        or throw("Can't open directory '$dir'");

	my @files = grep { /^[^.]/ && /${exp_regex}/ } readdir DIR;
	closedir DIR;

	die("No files match '$exp_regex' in input dir:\t$dir'") if (! @files);
	print "File regex '${exp_regex}' identified following input files:\n";

	system("mkdir -p $workdir/infiles") if (! -d "$workdir/infiles");

    foreach my $f (sort @files) {

	  next if (! -f "$dir/$f");
	  print "$f\t";
	  my $format = &get_file_format("$dir/$f");  #Maybe we can use some of the web code for this?
	  throw("File '$dir/$f' format is not a reconized format (e.g. bed|sam)") if ! $format;	  
	  print "\t$format file\n";

	  ### Check that files are gzipped
	  if(! &is_gzipped("$dir/$f")){
		
		if($zip){
		  my $cmd = "gzip $dir/$f";
		  print "Zipping input:\t$cmd\n";
		  warn "These gzip jobs needs submitting to the farm and waiting in the env CreateInputIDs";
 		  system($cmd) && die("Could not zip input:\t$dir/$f");
		  $f .='.gz';
		}
		else{
		  die("Input is not compressed with gzip:\t$dir/$f");
		}
	  }

	  
	  #why is this not matching CD4_DNase1.bed?
	  ###This always expects something like this ctype_ftype.SOMETHING.format.gz
	  #The SOMETHING is currently not being integrated into the exp name!
	  (my $experiment_name = $f) =~ s,(.*/)?([^/_]+_[^/_]+)\..*\.$format\.gz,$2,;
	  

		#This sets the last two tokens of the file name to the experiment name
	  #Should take this from the dir name? Or just set this manually 
	  #Could be based on the dir AUTHOR_PMID, but do this in the env and not here, so we can override
		#this with a defined experiment name.
		#And leave 
        $experiment_name .= '_'.$exp_suffix if ($exp_suffix);

		### validate cell and feature type
		my ($cell_type, $feature_type) = split('_', $experiment_name);
		throw ("Cell type '$cell_type' doesn't exist in database! ".
			   "Edit and rerun run_import_type.pl.") 
			unless (defined $cta->fetch_by_name($cell_type));
		throw ("Feature type '$feature_type' doesn't exist in database! ".
			   "Edit and rerun run_import_type.pl.") 
			unless (defined $fta->fetch_by_name($feature_type));

        # write input_id to database
        my $input_id = sprintf "%s:%s", $experiment_name, $f;
        #warn($input_id);

        # need to generate links in a workdir infiles directory to know where the 
        # files are that will be processed

	  my $link = "$workdir/infiles/$input_id";

	  if( ! -e $link){
        system("ln -s $dir/$f  $workdir/infiles/$input_id") && throw("Can't link 'ln -s $dir/$f  $workdir/infiles/$input_id'");
	  }

	  warn "Need to test for present input_ids here";

	  foreach my $lname(@logic_names){
		my $analysis_id = &get_analysis_id($lname);

		my $sql = "insert into input_id_analysis (input_id,analysis_id,input_id_type)".
		  " values ('${input_id}',${analysis_id},'file');";
		$pdb->dbc->do($sql);
		
		throw("Couldn't store input_id '$input_id'. Most likely it has already been ".
			  "stored. Drop your input_ids with CleanInputIds and rerun CreateInputIds.") if ($@);
	  }


	}
  } 
else {
  throw("Need to specify a valid input_type paramters e.g -slice, -dir or -array)");

}



sub get_analysis_id () {
  my $logic_name = shift;

    # get analysis_id of submit_type
    my $sql = "select analysis_id from analysis where logic_name='Submit${logic_name}'";
    
    my $sth = $pdb->dbc->prepare($sql);
    $sth->execute;

    throw("No analysis_id stored for logic_name='Submit${logic_name}'")
	  unless my $analysis_id = $sth->fetchrow;

    return $analysis_id;

}

1;
