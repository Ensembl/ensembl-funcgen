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

ensembl-efg rollback_array.pl

=head1 SYNOPSIS

rollback_array.pl [options]

Options:

Mandatory
  -arrays|a        List of Array names (space separated)
  -chips|c         List of ArrayChip design names to restrict the delete to
  -mode|m          Rollback mode supported by Helper::rollback_ArrayChip, default if full rollback i.e. 'probe'
  -pass|p          DB password
  -dbname|n        DB namme
  -port            DB port
  -host|h          DB host
  -user|u          DB user name
  -species         Latin name e.g. homo_sapiens
  -force|f         Performs a full delete, removing all records upto the specified level.  This is to specifically force the rollback of probe2transcript entries, as these can be associated with more than one array i.e. you may be deleting annotations which also belong to another array.
  -keep_xrefs      Flag to keep xrefs, even if we are deleting probe_features. Used when xrefs are external import i.e. not probe2transcript.
  #-no_clean_up     Skip optimize and and AUTO_INCREMENT reset. Useful if you want to save time when rolling 
                   back multiple sets of arrays.
  -help            Brief help message
  -man             Full documentation

=head1 DESCRIPTION

B<This program> removes all the imported probe, probe_feature and/or probe2transcript annotation data given set of Arrays or ArrayChips. This is not design to rollback arrays for which experimental data is present i.e. Nimblegen tiling arrays.

=cut



use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::Helper;
$| =1;

my (@array_names, $pass, $force, @chips, $vendor, %achips);
my ($host, $dbname, $species, $mode, $port, $user, $dnadb_pass, $keep_xrefs);
my ($dnadb_host, $dnadb_name, $dnadb_species, $dnadb_port, $dnadb_user);
my @tmp_args = @ARGV;


GetOptions (
			"arrays|a=s{,}"      => \@array_names,
			"chips|c=s{,}"    => \@chips,
			"vendor|v=s"      => \$vendor,
			"mode|m=s"        => \$mode,
			"dbuser|u=s"      => \$user,
			"dbpass|p=s"      => \$pass,
			"dbport=s"        => \$port,
			"dbname|n=s"      => \$dbname,
			"dbhost|h=s"      => \$host,
			"dnadb_pass=s"    => \$dnadb_pass,
			"dnadb_port=s"    => \$dnadb_port,
			"dnadb_name=s"    => \$dnadb_name,
			"dnadb_host=s"    => \$dnadb_host,
			"dnadb_user=s"    => \$dnadb_user,
			"species=s"       => \$species,
			"force|f"         => \$force,
			"keep_xrefs|k"    => \$keep_xrefs,
			"help|?"          => sub { pos2usage(-exitval => 0,  
												 -verbose => 2, 
												 -message => "Params are:\t@tmp_args");
									 },
		   ) or pod2usage(
						 -exitval => 1,
						 -message => "Params are:\t@tmp_args"
						);

die('Must define a -user parameter')    if ! $user;
die('Must define a -port parameter')    if ! $port;
die('Must define a -dbname parameter')  if ! $dbname;
die('Must define a -dbhost parameter')  if ! $host;
die('Must define a -pass parameter')    if ! $pass;
die('Must define a -species parameter') if ! $species;
$force = 'force' if $force;
$keep_xrefs = 'keep_xrefs' if $keep_xrefs;

my $dnadb;

if($dnadb_name){

  my $dnadb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
	(
	 -host => $dnadb_host || $host,
	 -dbname => $dnadb_name,
	 -user => $dnadb_user || $user,
	 -pass => $dnadb_pass || $pass,
	 -port => $dnadb_port,#don't default here as it may be on a different port
	 -species => $species,
	 -group => 'core'
	);
}

my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
													  -host => $host,
													  -dbname => $dbname,
													  -user => $user,
													  -pass => $pass,
						   						      -port => $port,
													  -species => $species,
													  -dnadb => $dnadb,
													  -group => 'funcgen',
													 );
$db->dbc->db_handle;#Test the DB connection



#This needs putting in Helper::rollback_Arrays
my $array_adaptor = $db->get_ArrayAdaptor;
my @arrays;


foreach my $aname(@array_names){
  my $array;

  $array = $array_adaptor->fetch_by_name_vendor($aname, $vendor);

  die ("Could not retrieve $vendor $aname Array") if ! $array;

  push @arrays, $array;
}

my $Helper = new Bio::EnsEMBL::Funcgen::Utils::Helper(
													  no_log => 1,#tees automatically with no_log
													 );


my %acs;

foreach my $array(@arrays){
  map {$acs{$_->design_id} = $_} @{$array->get_ArrayChips};
}

#do chips belong to experiment?
if(@chips){
  
  foreach my $chip_name(@chips){
	
	if(! exists $acs{$chip_name}){
	  die("$chip_name is not a valid ArrayChip design_id for the $vendor arrays:\t@array_names");
	}

	$achips{$chip_name} = $acs{$chip_name};	
  }
}
else{
  %achips = %acs;
}


$Helper->rollback_ArrayChips([values(%achips)], $mode, $force, $keep_xrefs);



1;
