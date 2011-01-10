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

create_input_ids.pl -- generates a list of input_ids from Encode regions for
    eFG analysis pipeline

=head1 SYNOPSIS

create_input_ids.pl OPTIONS

=head1 OPTIONS

    -host HOST        database host
    -port PORT        database port
    -user USER        user name
    -pass PASSWORD    password
    -dbname DBNAME    database name
    -slice            select slice as input id type
    -encode           uses encode regions as input_ids (w/ -slice)
    -toplevel         uses all toplevel slices as input_ids (w/ -slice)
    -array            select array as input id type
    -file             uses files in given directory (-dir) as input_ids
    -dir dir          directory to read infiles from (w/ -file)
    -exp_regex REGEX  regular expression to select certain experiments
                      (default: fetch all available if omitted)

=head1 DESCRIPTION

This script generates a list of input_ids from Encode regions / files for
setting up the eFG analysis pipeline.

=cut

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

my ($pass,$port,$host,$user,$dbname,$species,$data_version,
    $exp_regex,$exp_suffix,$slice,$encode,$toplevel,$array,$file,$dir,
	$help,$man,$debug, $submit_type);

$host = $ENV{EFG_HOST};
$port = $ENV{EFG_PORT};
$user = $ENV{EFG_WRITE_USER};
$dbname = $ENV{EFG_DBNAME};
$species = $ENV{SPECIES};
$data_version = $ENV{DATA_VERSION};

$exp_regex = '.*';


GetOptions (
            'pass|p:s'         => \$pass,
            'port:i'           => \$port,
            'host|h=s'         => \$host,
            'user|u=s'         => \$user,
            'dbname|d=s'       => \$dbname,
            'species=s'        => \$species,
            'data_version=s'   => \$data_version,
            'exp_regex|e=s'    => \$exp_regex,
            'exp_suffix=s'     => \$exp_suffix,
            'help|?'           => \$help,
            'man|m'            => \$man,
            'debug'            => \$debug,

            'slice'            => \$slice,
            'encode'           => \$encode,
            'toplevel'         => \$toplevel,

            'array'            => \$array,

            'file'             => \$file,
            'dir=s'            => \$dir,
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
throw("Must specify mandatory database hostname (-host).\n") if ! defined $host;
throw("Must specify mandatory database username. (-user)\n") if ! defined $user;
throw("Must specify mandatory database name (-dbname).\n") if ! defined $dbname;
throw("Must specify mandatory password (-pass).\n") if ! defined $pass;

$| = 1;

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info stack_trace_dump);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(open_file);
use Bio::EnsEMBL::Funcgen::Utils::Encode qw(get_encode_regions);

# Get eFG database adaptor

my $dnadb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (
     -host    => $ENV{CORE_HOST},
     -user    => $ENV{CORE_USER},
     -port    => $ENV{CORE_PORT},
     #-host    => 'ens-staging',
     #-user    => $user,
     #-port    => $port,
     -dbname  => $ENV{CORE_DBNAME},
     -species => $species,
     #-pass    => $pass,
     );

my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
    (
     -host    => $host,
     -user    => $user,
     -dbname  => $dbname,
#     -species => $species,
     -pass    => $pass,
     -port    => $port,
     -dnadb   => $dnadb,
     );

my $pdb = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
    (
     -host    => $host,
     -user    => $user,
     -dbname  => $ENV{'PDBNAME'},
     -species => $species,
     -pass    => $pass,
     -port    => $port,
     );

if ($slice) {
    
    $submit_type='Slice';

    # get analysis_id of submit_type
    my $analysis_id = &get_analysis_id($submit_type);
    
    # Get all experiments
    my $ea = $db->get_ExperimentAdaptor();
    my $exp = $ea->fetch_all();
        
    my @input_ids = ();

    if ($encode) {
        
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

    } elsif ($toplevel) {

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

        throw("Need to specify slice type (either -toplevel or -encode)");
        
    }

    foreach my $input_id (@input_ids) {
        my $sql = "insert into input_id_analysis (input_id,analysis_id,input_id_type)".
            " values ('${input_id}',${analysis_id},'${submit_type}');";
        #warn($sql);
        eval {
            $pdb->dbc->do($sql);
        };
        throw("Couldn't store input_id '$input_id'. Most likely it has already been ".
              "stored. Drop your input_ids with CleanInputIds and rerun CreateInputIds.")
            if ($@);
    }


} elsif ($array) {

    $submit_type='Array';

    # get analysis_id of submit_type
    my $analysis_id = &get_analysis_id($submit_type);

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
        my $sql = "insert into input_id_analysis (input_id,analysis_id,input_id_type)".
            " values ('${input_id}',${analysis_id},'${submit_type}');";
        #warn($sql);
        eval {
            $pdb->dbc->do($sql);
        };
        throw("Couldn't store input_id '$input_id'. Most likely it has already been ".
              "stored. Drop your input_ids with CleanInputIds and rerun CreateInputIds.")
            if ($@);
    }

} elsif ($file) {

    $submit_type='File';

	# get cell and feature type adapter
	my $cta = $db->get_CellTypeAdaptor();
	my $fta = $db->get_FeatureTypeAdaptor();

    # get analysis_id of submit_type
    my $analysis_id = &get_analysis_id($submit_type);
    
    if (! $dir) {
        throw("Need to specify a input directory containing ". 
              "files to be processed")
    }

    opendir(DIR, $dir)
        or throw("Can't open directory '$dir'");

    #$exp_regex='.*' unless ($exp_regex);

    my @files = grep { /^[^.]/ && /${exp_regex}/ } readdir DIR;

    closedir DIR;
  
    throw("No gzipped files found in input directory that match the regular expression '$exp_regex'")
        unless (@files);

    #print Dumper @files;

    unless (-d "$ENV{ANALYSIS_WORK_DIR}/infiles") {
        system("mkdir -p $ENV{ANALYSIS_WORK_DIR}/infiles");
    }

    foreach my $f (sort @files) {

        ### Check that files are gzipped
        throw("File is not compressed with gzip!") unless &is_gzip("$dir/$f");

        ### Check also that files are bed file format
        throw("File '$dir/$f' format is not bed format compliant!") unless &is_bed("$dir/$f");

        (my $experiment_name = $f) =~ s,(.*/)?([^/_]+_[^/_]+).*\.bed\.gz,$2,;
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
        system("ln -s $dir/$f  $ENV{ANALYSIS_WORK_DIR}/infiles/$input_id") == 0
            or throw("Can't link 'ln -s $dir/$f  $ENV{ANALYSIS_WORK_DIR}/infiles/$input_id'");

        my $sql = "insert into input_id_analysis (input_id,analysis_id,input_id_type)".
            " values ('${input_id}',${analysis_id},'${submit_type}');";
        #warn($sql);
        eval {
            $pdb->dbc->do($sql);
        };
        throw("Couldn't store input_id '$input_id'. Most likely it has already been ".
              "stored. Drop your input_ids with CleanInputIds and rerun CreateInputIds.")
            if ($@);

    } 

} else {

    throw("Need to specify type of input_ids to be specified ". 
          "(either -slice or -file)");
        
}

sub is_gzip {
    
    my ($file) = @_;

    open(FILE, "file -L $file |")
        or throw("Can't execute command 'file' on '$file'");
    my $retval = <FILE>;
    #print $retval, "\n";
    close FILE;

    return ($retval =~ m/gzip compressed data/) ? 1 : 0;
    
}

sub is_bed {
    
    my ($file) = @_;

    open(FILE, "zcat $file 2>&1 |")
        or throw("Can't execute command 'file' on '$file'");
    my @line = ();
    while (<FILE>) {
        chomp;
        @line = split("\t", $_);
        last;
    }
    close FILE;
    #print '<', join('><', @line),">\n";
    
    if (scalar @line < 6) {
        warn("Infile '$file' does not have 6 or more columns. We expect bed format: CHROM START END NAME SCORE STRAND.");
        return 0;
    #} elsif ($line[0] !~ m/^((chr)?[MTXYNT_\d]+)$/) {
    #    warn ("1st column must contain name of seq_region (e.g. chr1 or 1) in '$file'");
    #    return 0;
		#Commented this out for now due to HSCHR_RANDOM seqs
		#How does the webcode handle this?
    } elsif ($line[1] !~ m/^\d+$/ && $line[2] =~ m/^\d+$/) {
        warn ("2nd and 3rd column must contain start and end respectively in '$file'");
        return 0;
    } elsif ($line[5] !~ m/^[+-]$/) {
        warn ("6th column must define strand (either '+' or '-') in '$file'");
        return 0;
    }
    return 1;
    
}

sub get_analysis_id () {
    
    my ($type) = @_;

    # get analysis_id of submit_type
    my $sql = "select analysis_id from analysis where logic_name='Submit$type'";
    
    my $sth = $pdb->dbc->prepare($sql);
    $sth->execute;

    throw("No analysis_id stored for logic_name='Submit$type'")
        unless my $analysis_id = $sth->fetchrow;

    return $analysis_id;

}

1;
