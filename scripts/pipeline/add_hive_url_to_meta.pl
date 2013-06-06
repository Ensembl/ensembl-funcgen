#!/usr/bin/env perl

=head1 LICENSE


  Copyright (c) 1999-2013 The European Bioinformatics Institute and
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

add_hive_url_to_meta.pl - Sets/checks the hive url in the meta table

=head1 SYNOPSIS

add_hive_url_to_meta.pl PARAMETERS 

=head1 PARAMETERS

  Mandatory:
    -host     HOST        DB host
    -port     PORT        DB port (default = 3306)
    -user     USER        DB user name
    -pass     PASSWORD    DB password
    -dbname   DBNAME      DB name
    -url      HIVE_URL    The url of your hive DB e.g. 
                            'mysql://DB_USER:DB_PASS@DB_HOST:DB_PORT/DB_NAME'
    -pipeline PIPELINE    Name of the pipeline (i.e. $ENV_NAME from your derivative pipeline.env)

  Optional:
    -species  SPECIES     Only required where species is not already set in meta or for multi species DBs
    -help
    -man

=head1 DESCRIPTION

This script checks whether a hive url meta tabe entry has been set for the 
given pipeline. If it is absent it sets it, if it does not if dies, if it is 
absent, it sets the value.

The script is used during initialisation of the pipeline environment. The aim 
is to ensure that it is not possible to run two versions of the same pipeline
on the same output DB.

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;


my ($pass, $host, $user, $dbname, $species, $url, $pname, $help, $man);
my $port = 3306;
my @tmp_args = @_;

GetOptions (
            'pass:s'     => \$pass,
            'port:i'     => \$port,
            'host=s'     => \$host,
            'useru=s'    => \$user,
            'dbname=s'   => \$dbname,
            'species=s'  => \$species,
            'url=s'      => \$url,
            'pipeline=s' => \$pname,
      
            'help'             => sub { pos2usage(-exitval => 0); }, #do we need verbose here?
            #removed ~? frm here as we don't want to exit with 0 for ?
            
            'man|m'            => sub { pos2usage(-exitval => 0,  -verbose => 2); },
           ) or pod2usage(-exitval => 1, -message => "Specified parameters are:\t@tmp_args"); 
           


### check options ###
if( ! (defined $user &&
       defined $pass &&
       defined $dbname &&
       defined $host &&
       defined $url &&
       defined $pname) ){
  pod2usage(-exitval => 1, 
            -message => "You must specify all mandatory parameters\nSpecified parameters are:\t@tmp_args"
            );        
 }

throw("Must specify mandatory database hostname (-host).\n") if ! defined $host;
throw("Must specify mandatory database username. (-user)\n") if ! defined $user;
throw("Must specify mandatory database name (-dbname).\n") if ! defined $dbname;
throw("Must specify mandatory password (-pass).\n") if ! defined $pass;

$| = 1;


#This will currently fail if we don't have a dnadb on live-mirror
#todo add a -no_dnadb to the DBADaptor

my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
    (
     -host    => $host,
     -user    => $user,
     -dbname  => $dbname,
     -pass    => $pass,
     -port    => $port,
     -species => $species,
 #    -dnadb   => $dnadb,
     );

#Test the connection 
$db->dbc->db_handle;

my $mc         = $db->get_MetaContainer;
my $meta_key   = $pname.'_hive_url';
my $meta_value = $mc->single_value_by_key($meta_key);

if( ! defined $meta_value){ #Store the new meta entry
  #Need to tie this to species_id, such that we don't get NULL in unique key
  #Just need to validate we have species set
  #add multi species support  
  
  my $db_species = $db->species;
  
  #sanity check
  if(! defined $db_species){
    die("No species defined in the DB. Unsafe to add meta entry due to NULL in unique key".
      "\nPlease add manually");
  }
  elsif(defined $species &&
        ($species ne $db_species) ){
    die("The species defined in the DB($db_species) does not match the species parameter ($species)");
  }
  
  #Add key via API to store with appropriate species_id
  eval { $mc->store_key_value($meta_key, $url); };
  
  if($@){
    die("Failed to store hive url meta entry.\n$@");  
  }

}
elsif($meta_value ne $url){
  die("$dbname is currently locked to a different $pname hive DB:\t$meta_value\n".
    'Please use that hive DB or drop that pipeline '.
    'i.e. DropPipeline in that pipeline environment or remove meta entry');
}
#else all is good.

1;
