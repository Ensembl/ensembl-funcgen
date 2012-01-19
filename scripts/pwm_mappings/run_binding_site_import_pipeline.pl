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

run_binding_site_import_pipeline.pl -- deletes previous data and runs the binding site import

=head1 DESCRIPTION

Deletes all motif features for matrices that we're importing and sets up the import pipeline for import

=head1 OPTIONS

=over

=item B<help>

Give short help

=item B<-dbhost>

Host where the database is

=item B<-dbuser>

User of the database 

=item B<-dbpass>

Password for the database user 
Needs writing access

=item B<-dbport>

Port of the host where the database is 

=item B<-dbname>

Name of the database 


=item B<-dnadbhost>

Host of the specific core database to use 

=item B<-dnadbuser>

User of the specific core database 

=item B<-dnadbport>

Port of the host where the specific core database to use is

=item B<-dnadbname>

Name of the specific core database to use

=item B<-workdir>

Folder where the data is found

=item B<-output_dir>

Folder to output the results of the pipeline

=item B<-slices>

List of slices to be imported, separated by ; eg. 1;2;X;23

=item B<-continue>

If set, it will assume a pipeline db already exists

=back

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

my ($host, $port, $user, $pass, $dbname);
my ($dnadb_host, $dnadb_port, $dnadb_user, $dnadb_name);
my ($help, $workdir, $output_dir, @slices);

#get command line options
print "run_binding_site_import.pl @ARGV\n";

GetOptions (
	    'dnadb_host=s'       => \$dnadb_host,
	    'dnadb_user=s'       => \$dnadb_user,
	    'dnadb_port=i'       => \$dnadb_port,
	    'dnadb_name=s'       => \$dnadb_name,
	    'dbhost=s'           => \$host,
	    'dbuser=s'           => \$user,
	    'dbport=i'           => \$port,
	    'dbpass=s'           => \$pass,
	    'dbname=s'           => \$dbname,
	    'workdir=s'          => \$workdir,
	    'output_dir=s'       => \$output_dir,
	    'slices=s{,}'        => \@slices,
	    "help|h"             => \$help,
	   )  or pod2usage( -exitval => 1 ); #Catch unknown opts

pod2usage(1) if ($help);

if(!$host || !$port || !$user || !$dbname ) {  print "Missing connection parameters for efg db\n"; pod2usage(0); }
if(!$dnadb_host || !$dnadb_port || !$dnadb_user || !$dnadb_name ) {  print "Missing connection parameters for core db\n"; pod2usage(0); }
if(!$workdir) {  print "Missing working folder(s)\n"; pod2usage(0); }
if(!$output_dir) {  print "Missing output folder(s)\n"; pod2usage(0); }

my $coredb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
  (
   '-host'        => $dnadb_host,
   '-user'        => $dnadb_user,
   '-port'        => $dnadb_port,
   '-dbname'      => $dnadb_name,
  );

my $efgdb = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
  (
   -host    => $host,
   -port    => $port,
   -user    => $user,
   -dbname  => $dbname,
   -pass    => $pass,
   -dnadb  => $coredb,
  );

# Test connection
$efgdb->dbc->db_handle;

opendir(DIR,$workdir);
my @files = readdir(DIR);
closedir DIR;

my $first = 1;

my $bma = $efgdb->get_BindingMatrixAdaptor();

foreach my $file (@files){
  
  next if(!( $file =~ /^(.*).pwm_map$/));
  my $matrix = $1;
  print "Matrix: ".$matrix."\n";

  my @bms = @{ $bma->fetch_all_by_name($matrix) };
  if(scalar(@bms)==0) { warn $matrix." does not exist - ignoring"; next; }
  if(scalar(@bms)>1) { warn "More than onw matrix with name ".$matrix." - ignoring"; next; }

  my $bm = $bms[0];

  my $mf_sql =  "motif_feature where binding_matrix_id=".$bm->dbID;
  if(scalar(@slices)>0){
    warn "Restricting to requested ".scalar(@slices)." slices";
    $mf_sql = $mf_sql." and seq_region_id in (select seq_region_id from seq_region where name in ('".join("','",@slices)."'))";
  }

  warn "Deleting previous data for $matrix with matrix id ".$bm->dbID." - this will likely imply that a new regulatory build is needed";
  $efgdb->dbc->do("delete from associated_motif_feature where motif_feature_id in (select motif_feature_id from ".$mf_sql." )");
  $efgdb->dbc->do("delete from ".$mf_sql);

  my $cmd="init_pipeline.pl Bio::EnsEMBL::Funcgen::HiveConfig::ImportMotifFeatures_conf -dnadb_host $dnadb_host -dnadb_port $dnadb_port -dnadb_user $dnadb_user -dnadb_name $dnadb_name -host $host -port $port -user $user -pass $pass -dbname $dbname -output_dir $output_dir -efg_src $ENV{SRC}/ensembl-functgenomics/ -file ${workdir}/${file} -matrix $matrix ".(scalar(@slices)>0 ? " -slices ".join(",",@slices) : "")." ".($first ? '' : " -job_topup");
  print $cmd."\n";
  system($cmd);

  $first = 0;

}


