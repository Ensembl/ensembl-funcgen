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

generate_new_Stable_ids.pl

=head1 SYNOPSIS

This script generates new stable IDs for the first build of a given species.

=head1 OPTIONS

 Mandatory:
  -species      Latin name as used in DB name or meta table e.g. homo_sapiens
  -host         efg DB host(assumes core host is ens-staging)
  -schema_build The schema_build of the DB you want to update e.g. 54_37g
  -pass         The DB password for the efg DB

 Optional:
  -help         Prints this POD and exits.


=cut


BEGIN{
	if(! defined $ENV{'EFG_DATA'}){
		if(-f "$ENV{'HOME'}/src/ensembl-functgenomics/scripts/.efg"){
			system (". ~/src/ensembl-functgenomics/scripts/.efg");
		}else{
			die ("This script requires ensembl-functgenomics/scripts/.efg\n".
				 "Please source it before running this script\n");
		}
	}
}

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Bio::EnsEMBL::Funcgen::Utils::Helper;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(species_chr_num open_file);


$| = 1;
my ($species, $schema_build, $pass, $dbhost);
my @tmp_args = @ARGV;

GetOptions(
		   'species=s'      => \$species,
		   'schema_build=s' => \$schema_build,
		   'pass=s'         => \$pass,
		   'host=s'         => \$dbhost,
		   'help'                   => sub { pos2usage(-exitval => 0, -message => "Params are:\t@tmp_args"); }
	   ) or pod2usage(
					  -exitval => 1,
					  -message => "Params are:\t@tmp_args"
					 );


if(! ($species && $schema_build && $pass && $dbhost)){
  pod2usage(
			-exitval => 1,
			-message => "You have failed to specifiy a mandatory parameter\nParams are:\t@tmp_args"
		   );
  exit 1;
}
#Let's just hard code this for staging as we will always be present after the build


my $cdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
  (
   -host => 'ens-staging',
   -dbname => $species.'_core_'.$schema_build,
   -species => $species,
   -user => "ensro",
  );
$cdb->dbc->db_handle;#Test

my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
  (
   -host => $dbhost,
   -dbname => $species.'_funcgen_'.$schema_build,
   -species => $species,
   -dnadb => $cdb,
   -user => 'ensadmin',
   -pass => $pass,
  );
$db->dbc->db_handle;#Test



my $fset = $db->get_FeatureSetAdaptor->fetch_by_name('RegulatoryFeatures');

#Test we have a FeatureSet and no previously assigned stable IDs
if(! $fset){
  die("No RegulatoryFeatures FeatureSet available in ${species}_funcgen_${schema_build}");
}


my $cmd = 'SELECT stable_id from regulatory_feature rf where feature_set_id='.$fset->dbID.' order by stable_id desc limit 1';
my ($next_stable_id) = @{$db->dbc->db_handle->selectrow_arrayref($cmd)};

if($next_stable_id){
  die('There are already stable IDs present in the regulatory_feature table, maybe you want to run stable_id_mapper.pl?');
}


#Now assign new IDs
#To ensure we utilise all available IDs we need to do incremental update
$cmd = 'SET @count=0';
$db->dbc->do($cmd);
$cmd = 'update regulatory_feature set stable_id=(SELECT @count:=@count + 1) where feature_set_id='.$fset->dbID;
my $count = $db->dbc->do($cmd);

print "Updated $count stable IDs\n";



1;

