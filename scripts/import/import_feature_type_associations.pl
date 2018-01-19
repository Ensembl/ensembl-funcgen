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

ensembl-efg import_feature_type_associations.pl
  
=head1 SYNOPSIS

import_feature_type_associations.pl [options]

Options:

Mandatory

  -species|s       Species name

  -pass|p        The password for the EFG DB

  -dbname          Defines the eFG dbname

  -file            File containing the associations
  
Optional

  -host|h        Defines the eFG db host [ens-genomics1]

  -port|l        Defines the eFG db port [3306]

  -user|u        Defines the eFG db user [ensadmin]

  -help            Brief help message

  -man             Full documentation

=head1 OPTIONS

=over 8

=item B<-species|s>

Species name

=item B<-pass|p>

Database password

=item B<-dbname>

Database name

=item B<-host|h>

Database host (defaults to ens-genomics1)

=item B<-port|l>

Database port (defaults to 3306)

=item B<-user|u>

Database user (defaults to 'ensadmin')

=item B<-dnadb_host>

Core Database host (defaults to ens-livemirror)

=item B<-dnadb_port>

Core Database port (defaults to 3306)

=item B<-dnadb_user>

Core Database user (defaults to 'ensro')

=item B<-file>

Input file with data

=item B<-help>

Print a brief help message and exits.

=back

=head1 DESCRIPTION

B<This script> adds associations between feature types.
This is tipically to link Transcription Factor Complexes to their components

=cut

use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );

my ($dnadb_host,$dnadb_port,$dnadb_user);
my ($pass,$port,$host,$user,$dbname,$species,$file,$help);
$port = 3306;
$host='ens-genomics1';
$user='ensadmin';

GetOptions (
	    "dnadb_host=s" => \$dnadb_host,
	    "dnadb_port=s" => \$dnadb_port,
	    "dnadb_user=s" => \$dnadb_user,
	    "pass|p=s"     => \$pass,
	    "port|l=s"     => \$port,
	    "host|h=s"     => \$host,
	    "user|u=s"     => \$user,
	    "dbname=s"     => \$dbname,
	    "species|s=s"  => \$species, 
	    "file|f=s"     => \$file, 
	    "help|?"       => \$help,
	   );
pod2usage(1) if $help;
pod2usage(0) if !defined($species);
pod2usage(0) if !defined($dbname);
pod2usage(0) if !defined($pass);
pod2usage(0) if !defined($file);

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
$efg_dba->dbc->db_handle;

print ":: Loading FeatureType Associations from file $file\n";
my $fta	= $efg_dba->get_FeatureTypeAdaptor();     

my %fts;
open(FILE,$file);
<FILE>; #title
while(<FILE>){
	chomp;
	my ($feature_type, $assoc_ft) = split(/\s+/);
	my $ft1 = $fta->fetch_by_name($feature_type);
	if(!$ft1){ warn $feature_type." does not exist"; next; }
	my $ft2 = $fta->fetch_by_name($assoc_ft);
	if(!$ft2){ warn $assoc_ft." does not exist"; next; }
	push @{$fts{$feature_type}}, $assoc_ft;
}
close FILE;

foreach my $ft (keys %fts){
  #print $ft."\n";
  my $ft_obj = $fta->fetch_by_name($ft);
  #print $ft1->name."\n";
  my @afts = map { $fta->fetch_by_name($_); } @{$fts{$ft}};
  #map { print "\t".$_->name."\n"; } @afts;
  eval{
    $ft_obj->associated_feature_types(\@afts);
    $fta->store_associated_feature_types($ft_obj);
    #Do the reverse association...
    map{
      $_->associated_feature_types([ $ft_obj ]);
      $fta->store_associated_feature_types($_);
    } @afts;
  };
  if($@){ warn $@; };
}

