#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=head1 NAME

update_DB_for_release.pl
  

=head1 SYNOPSIS

This script performs several updates to the eFG DBs as part of the release cycle.

 perl update_DB_for_release.pl <mandatory paramters> [ optional parameters ]


=head1 DESCRIPTION

This script uses the Bio::EnsEMBL::Funcgen::Utils::Healthchecker module to perform critical updates and checks on the funcgen DBs during the release cycle. 


This default behaviour can be over-ridden by specifying a space separated list of method names e.g.

  -methods validate_new_seq_regions update_meta_schema_version check_meta_strings analyse_and_optimise_tables set_current_coord_system update_meta_coord clean_xrefs log_data_sets check_stable_ids


See Bio::EnsEMBL::Funcgen::Utils::HealthChecker for more details.


=head1 OPTIONS

 Mandatory
  -species      Latin name as used in DB name or meta table e.g. homo_sapiens
  -host         DB host
  -user         DB user
  -pass         DB pass
  -data_version e.g. 58_37k
  

 Optional
  -dbname            DB name (default production name will be used if not specified)
  -dnadb_host        DNADB host
  -skip_meta_coord   These skip methods are useful as these are likely to compete
  -skip_analyse      with running processes, but you may want to run the rest of the
  -skip_xref_cleanup update, saving these skipped methods until everything is finished.
  -check_displayable Forces log_data_sets to only use DISPLAYABLE sets
  -methods           Over-rides defaults behaviour and only run specified methods (see above)
  -tee               Tees output to STDOUT
  -log_file           
  -no_log            No log file, but turns on tee
  -help              Prints a short help page
  -man               Prints the entire POD documentation


=head1 EXAMPLE

perl update_DB_for_release.pl -species homo_sapiens -host ens-genomics1 -user $USER -data_version 58_37c -dbname homo_sapiens_funcgen_58_37c -dnadb_host ens-staging -check_displayable -pass $PASS

=head1 SEE ALSO

Bio:EnsEMBL::Funcgen::Utils::HealthChecker

=cut



use strict;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::HealthChecker;
use Pod::Usage;
use Getopt::Long;

my ($pass, $species, $schema_build, $skip_meta_coord, $dnadb_host, $dnadb_user, $dnadb_pass);
my ($dnadb_port, $check_displayable);
my ($help, $man, $dbname, $skip_xref_cleanup, $skip_analyse, @methods, @method_params);
my $user = 'ensadmin';
my $port = 3306;
my $host = 'ens-genomics1';
my @tmp_args = @ARGV;

#define here to avoid warnings
$main::_tee = 0;

GetOptions( 
		   'host|h=s'          => \$host,
		   'port=i'            => \$port,
		   'user|u=s'          => \$user,
		   'pass|p=s'          => \$pass,
		   'dbname=s'          => \$dbname,
		   'species|d=s'       => \$species,
		   'data_versions=s'   => \$schema_build,#mandatory as default ensembldb will not exist
		   'dnadb_host=s'      => \$dnadb_host,
		   'dnadb_user=s'      => \$dnadb_user,
		   'dnadb_pass=s'      => \$dnadb_pass,
		   'dnadb_port=s'      => \$dnadb_port,
		   'skip_meta_coord'   => \$skip_meta_coord,
		   'skip_analyse'      => \$skip_analyse,
           'skip_xref_cleanup' => \$skip_xref_cleanup,
		   'check_displayable' => \$check_displayable,
		   'methods=s{,}'      => \@methods,
		   'method_params=s{,}'=> \@method_params,
		   'help|?'            => \$help,
		   'man|m'             => \$man,
		   'log_file=s'        => \$main::_log_file,
		   'no_log'            => \$main::_no_log,
		   'tee'               => \$main::_tee,
		   #'slice=s'          => \$test_slice,
		   #skip dumps?
		   #force update
		   #add opt for old, new & stable fset name
		  ) or pod2usage(
						 -exitval => 1,
						 -message => "Params are:\t@tmp_args"
						);

pod2usage(0) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;


my @builds = @ARGV;

$dbname = "${species}_funcgen_${schema_build}" if ! defined $dbname;

if(! $main::_no_log){
  $main::_log_file ||= $ENV{'HOME'}."/logs/update_DB_for_release.$dbname.$$.log";
  print "Writing log to:\t".$main::_log_file."\n";
}


my $efg_db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
														  -dbname  => $dbname,
														  -host    => $host,
														  -port    => $port,
														  -user    => $user,
														  -pass    => $pass,
														  -species => $species,
														  -dnadb_host => $dnadb_host,
														  -dnadb_user => $dnadb_user,
														  -dnadb_pass => $dnadb_pass,
														  -dnadb_port => $dnadb_port,
														 );

#Test the db connections
$efg_db->dbc->db_handle;
$efg_db->dnadb->dbc->db_handle;

my $hchecker = Bio::EnsEMBL::Funcgen::Utils::HealthChecker->new(
																-db                => $efg_db,
																-builds            => \@builds,
																-skip_meta_coord   => $skip_meta_coord,
																-skip_analyse      => $skip_analyse,
																-skip_xref_cleanup => $skip_xref_cleanup,
																-check_displayable => $check_displayable,
															   );






if(@methods){
  foreach my $method(@methods){

	if(! $hchecker->can($method)){
	  die("You have passed an invalid method:t\$method");
	}

	$hchecker->$method;  
  }
}
elsif(@method_params){
  my ($method, @params);
  ($method, @params) = @method_params;

  if(! $hchecker->can($method)){
	die("You have passed an invalid method:t\$method");
  }

  $hchecker->$method(@params); 
}
else{
 $hchecker->update_db_for_release();

}

