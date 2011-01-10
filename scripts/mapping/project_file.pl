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

project_file.pl - Projects a features in flat file to another assembly

=head1 SYNOPSIS

project_file.pl [arguments]

Required arguments:

  --dbname, db_name=NAME              database name NAME
  --host, --dbhost, --db_host=HOST    database host HOST
  --port, --dbport, --db_port=PORT    database port PORT
  --user, --dbuser, --db_user=USER    database username USER
  --pass, --dbpass, --db_pass=PASS    database passwort PASS

Optional arguments:

  --conffile, --conf=FILE             read parameters from FILE
                                      (default: conf/Conversion.ini)

  --logfile, --log=FILE               log to FILE (default: *STDOUT)
  --logpath=PATH                      write logfile to PATH (default: .)
  --logappend, --log_append           append to logfile (default: truncate)
  --loglevel=LEVEL                    define log level (default: INFO)

  --is_component, --is-component      script is called from a wrapper script

  -i, --interactive                   run script interactively (default: true)
  -n, --dry_run, --dry                don't write results to database
  -h, --help, -?                      print help (this message)

=head1 DESCRIPTION

This script projects features from an old assembly to a new assembly

=cut

use strict;
use warnings;
#no warnings 'uninitialized';

use FindBin qw($Bin);
use Bio::EnsEMBL::Utils::ConfParser;
use Bio::EnsEMBL::Utils::Logger;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw (open_file);

use Bio::EnsEMBL::Funcgen::Utils::Helper;#replace logger or inherit from logger?

#This is quite useful, lists params defined
#Summarises num wanring runtime, complete time and memusage at end.
# parse configuration and commandline arguments
my $conf = new Bio::EnsEMBL::Utils::ConfParser(
  -SERVERROOT => "$Bin/../../..",
  -DEFAULT_CONF => ""
);

$conf->parse_options
  (
   'host=s' => 1,
   'port=n' => 0,
   'user=s' => 1,
   'pass=s' => 0,
   'dbname=s' => 1,
   'species=s' => 0,
   'ignore_length' => 0,
   #'force_store' => 1,
   'file=s' => 1,
   'old_assembly=s' => 1,
   'new_assembly=s' => 1,
   'coord_system=s' => 0,
   'associations'      => 0,
  );


#Should use Bed parser here!

$main::_no_log = 1;
#$main::_tee    = 1;#Isn't this set by default if no log i set?
my $helper = new Bio::EnsEMBL::Funcgen::Utils::Helper;

#assign assemblies for regex embedding
my $new_assembly = $conf->param('new_assembly');
my $old_assembly = $conf->param('old_assembly');
#my $cs_level     = $conf->param('coord_system') || 'chromosome';


# initialise log


# connect to database and get adaptors

#We only need the old DB to get access to the old toplevel
#So can use ensembldb
#This may include contigs which have now been merged into a chromosome
#and so would not be easily accessable via the new DB
#we would have to somehow get all non toplevel regions which are not 
#included in the new assembly.
#But they would not be present in the new DB
#We don't actually have this cross level mapping anyway, so can remove
#And only use new DB!

#Need to expose old_db params here!

my $efg_db = new Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor
  (
   -host   => $conf->param('host'),
   -port   => $conf->param('port'),
   -user   => $conf->param('user'),
   -pass   => $conf->param('pass'),
   -dbname => $conf->param('dbname'),
   -species =>  $conf->param('species'),
   -group  => 'funcgen',
   -dnadb_host => 'ens-staging',
   -dnadb_user => 'ensro',
   #-dnadb  => $new_cdb,
  );



#Test connections
$efg_db->dbc->db_handle;
$efg_db->dnadb->dbc->db_handle;
my $sa = $efg_db->dnadb->get_SliceAdaptor;
my $file = $conf->param('file');

my $fh   = open_file($file);
my %strands = (
			  # '+' => 1,
			  # '-' => -1,
			   'F' => 1,
			   'R' => -1,
			   
			   
			   '.' => 0,
			  );



my $outfile = $file.'.out';
$outfile = open_file($outfile, '>');
my %counts;

foreach my $line(<$fh>){

  #bed format
  chomp $line;
  #my ($chr, $start, $end, $name, $score, $strand, @other_fields) = split/\s+/o, $line;
  my ($chr, $start, $strand, $seq) = split/\s+/o, $line;
  
  my $end = ($start + length($seq) -1 );
  my $score = '';
  my $name = '';
  my @other_fields;
  $strand = $strands{$strand};
  

  $chr =~ s/chr//;
  
  $chr ='MT' if ($chr eq 'M');
  $start +=1; #UCSC half open!

  my $slice = $sa->fetch_by_region(undef, $chr, $start, $end, $strand, $old_assembly);
  #project using =ve strand to simplyfy local sqe_start/end calcs

  if(! $slice){
	warn "Could not get slice $old_assembly:$chr:$start:$end:$strand";
	next;
  }


  my @segments = @{$slice->project($slice->coord_system->name, $new_assembly)};

  if(scalar(@segments) == 1){
	my ($seg_start, $seg_end, $seg_slice) = @{$segments[0]};

	#seq_start and end are local to slice!!!
	#We need to take account of strand issues here?
	#Let's just project the +ve strand to simplify the calcs
	#And then just sub the srtand back in?
	#Unless our returned slice is not on the same strand!
	if($seg_slice->strand != 1){
	  $start = ($seg_slice->end - $seg_end + 1);
	  $end   = ($seg_slice->end - $seg_start +1);
	  #123456
	  #  4  1
	}
	else{
	  $start = ($seg_start + $seg_slice->start -1);
	  $end   = ($seg_end + $seg_slice->start - 1);
	}

	$counts{$chr}{mapped}++;

	print $outfile join("\t", ($seg_slice->seq_region_name, ($seg_start + $seg_slice->start -1), ($seg_end + $seg_slice->start - 1), $name, $score, $strand, @other_fields))."\n";
  }
  else{
	warn "Failed to project:\t$line\n";
	$counts{$chr}{fail}++;
  }
}

my $mapped = 0;
my $failed = 0;

foreach my $chr(keys %counts){
  print "Projected features on chromsome $chr:\t".$counts{$chr}{mapped}."\n";
  print "Failed feature projections on chromsome $chr:\t".($counts{$chr}{fail} || 0)."\n";
  $mapped += $counts{$chr}{mapped};
  $failed += $counts{$chr}{fail};
}


print "Project $mapped features from $old_assembly to $new_assembly
Failed to map $failed feature\n";

#Add sort here
#Use parser sort key? Or EFGUtils

close $file;
close $outfile;
