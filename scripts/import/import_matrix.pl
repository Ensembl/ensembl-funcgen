#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

import_matrix.pl
  
=head1 SYNOPSIS

import_matrix.pl [options]

Options:

Mandatory

  -species|s       Species name

  -matrix_id|m     Jaspar Id of the Matrix

  -feature_type|f       Feature Type name

  -matrix_data_dir|d     Folder where Jaspar data (matrix_only.txt) is

  -dbpass|p        The password for the EFG DB

  -dbname          Defines the eFG dbname
  
Optional

  -dbhost|h        Defines the eFG db host [ens-genomics1]

  -dbport|l        Defines the eFG db port [3306]

  -dbuser|u        Defines the eFG db user [ensadmin]

  -help            Brief help message

  -man             Full documentation

=head1 OPTIONS

=over 8

=item B<-species|s>

Species name

=item B<-matrix_data_dir>

Folder where the Jaspar matrix data is (matrix_only.txt)

=item B<-matrix_id|m>

Jaspar Id of the matrix to be uploaded (e.g. 'MA0139.1')

=item B<-feature_type|f>

Feature Type to be associated to the matrix (e.g. 'CTCF')
This feature type needs to exist in the database

=item B<-dbpass|p>

Database password

=item B<-dbname>

Database name

=item B<-dbhost|h>

Database host (defaults to ens-genomics1)

=item B<-dbport|l>

Database port (defaults to 3306)

=item B<-dbuser|u>

Database user (defaults to 'ensadmin')

=item B<-dnadb_host>

Core Database host

=item B<-dnadb_port>

Core Database port

=item B<-dnadb_user>

Core Database user

=item B<-help>

Print a brief help message and exits.

=back

=head1 DESCRIPTION

B<This script> Imports a specific matrix and associates it to a specific Feature Type

=cut

use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Funcgen::Utils::Helper;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::BindingMatrix;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );

my ($dnadb_host,$dnadb_port,$dnadb_user, $dnadb_name);
my ($pass,$port,$host,$user,$dbname,$species,$matrix_data_dir,$matrix_id,$feature_type,$help);
$port = 3306;
$host='ens-genomics1';
$user='ensadmin';

GetOptions (
	    "dnadb_host=s" => \$dnadb_host,
	    "dnadb_port=s" => \$dnadb_port,
	    "dnadb_user=s" => \$dnadb_user,
	    "dnadb_name=s" => \$dnadb_name,
	    "dbpass|p=s"    => \$pass,
	    "dbport|l=s"    => \$port,
	    "dbhost|h=s"    => \$host,
	    "dbuser|u=s"    => \$user,
	    "dbname=s"      => \$dbname,
	    "species|s=s"   => \$species, 
	    "matrix_id|m=s" => \$matrix_id,
	    "matrix_data_dir|d=s" => \$matrix_data_dir,
	    "feature_type|f=s" => \$feature_type,
	    "help|?"        => \$help,
	   );
pod2usage(1) if $help;
pod2usage(0) if !defined($species);
pod2usage(0) if !defined($matrix_id);
pod2usage(0) if !defined($feature_type);
pod2usage(0) if !defined($dbname);
pod2usage(0) if !defined($pass);

my $efg_dba = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
    (
     -host    => $host,
     -user    => $user,
     -dbname  => $dbname,
     -species => $species,  
     -port    => $port,
     -pass    => $pass,
     -dnadb_host => $dnadb_host,
     -dnadb_port => $dnadb_port,
     -dnadb_user => $dnadb_user,
     );


my $aa	= $efg_dba->get_AnalysisAdaptor();
my $analysis = $aa->fetch_by_logic_name('Jaspar');
if(!$analysis){ throw "An Analysis named Jaspar needs to be in the database"; }

my $fta	= $efg_dba->get_FeatureTypeAdaptor();     
my $ft = $fta->fetch_by_name($feature_type);
if(!$ft){ throw "Could not find Feature Type $feature_type"; }

my $bma	= $efg_dba->get_BindingMatrixAdaptor();

open(FILE,$matrix_data_dir."/matrix_only.txt");
while(<FILE>){
  if(/^>$matrix_id/){

    #my $matrix = $_.<FILE>.<FILE>.<FILE>.<FILE>;
    #ignore the header from now on...
    my $matrix = <FILE>.<FILE>.<FILE>.<FILE>;

    if(@{$bma->fetch_all_by_name_FeatureType($matrix_id,$ft)}){ 
      warn "Matrix already exists in the database. Skipping import"; 
    } else {
      my $bm = Bio::EnsEMBL::Funcgen::BindingMatrix->new(
							 -name         => $matrix_id,
							 -analysis     => $analysis,
							 -feature_type => $ft,
							 -description  => 'Jaspar Matrix',
							 -frequencies  => $matrix
							);
      print $bm->frequencies."\n";
      print $bm->frequencies_revcomp."\n";

      ($bm) = $bma->store($bm);
      if(!$bm){ 
	warn "Failure on import"; 
	close FILE;
	exit 1;
      } else {
	warn "Matrix imported successfully"; 
	close FILE;
	exit 0;
      }
    }
  }
}
close FILE;

warn "Could not find $matrix_id";
exit 1;
