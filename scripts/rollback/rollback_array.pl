#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

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
use Data::Dumper;
$| =1;

my (@array_names, $pass, $force, @chips, $vendor);
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
    "help|?"          => sub { 
      pos2usage(
        -exitval => 0,  
        -verbose => 2, 
        -message => "Params are:\t@tmp_args"
      );
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
$db->dbc->db_handle || die("Can't connect to funcgen database!");

my $array_adaptor = $db->get_ArrayAdaptor;
my @arrays;
foreach my $array_name (@array_names){

  my $array = $array_adaptor->fetch_by_name_vendor($array_name, $vendor);
  if (!$array) {
      warn ("Could not retrieve $vendor $array_name Array");
      next;
  }
  push @arrays, $array;
}

my $Helper = new Bio::EnsEMBL::Funcgen::Utils::Helper(
  no_log => 1,
);

my %array_chip_designId_to_arrayChip;

foreach my $array (@arrays){
  map {$array_chip_designId_to_arrayChip { $_->design_id } = $_ } @{$array->get_ArrayChips};
}

my %array_chips;
if(@chips){  
  foreach my $chip_name(@chips){
	
	if(! exists $array_chip_designId_to_arrayChip{$chip_name}){
	  die("$chip_name is not a valid ArrayChip design_id for the $vendor arrays:\t@array_names");
	}
	$array_chips{$chip_name} = $array_chip_designId_to_arrayChip{$chip_name};	
  }
} else{
  # Design ids in array_chip seem to be the same as the name in the array table.
  %array_chips = %array_chip_designId_to_arrayChip;
}

$Helper->rollback_ArrayChips([values(%array_chips)], $mode, $force, $keep_xrefs);

1;
