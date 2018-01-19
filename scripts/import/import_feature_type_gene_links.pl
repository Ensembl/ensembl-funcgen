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

ensembl-efg import_feature_type_gene_links.pl
  
=head1 SYNOPSIS

import_feature_type_gene_links.pl [options]

Options:

Mandatory

  -species|s       Species name

  -pass|p        The password for the EFG DB

  -dbname          Defines the eFG dbname

  -file            Specifies file with data
  
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

Core Database host

=item B<-dnadb_port>

Core Database port

=item B<-dnadb_user>

Core Database user


=item B<-file|f>

Data file (usually in types/$species.FeatureType_Genes.txt)

=item B<-help>

Print a brief help message and exits.

=back

=head1 DESCRIPTION

B<This script> adds ENSEMBL Gene Ids as external references to Feature Types.
It reads the information from a text file with genes linked to the Feature Types
This file is located in types/$species.FeatureType_Genes.txt

=cut

use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::Funcgen::Utils::Helper;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(add_external_db);
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );

my ($dnadb_host,$dnadb_port,$dnadb_user, $dnadb_name);
my ($pass,$port,$host,$user,$dbname,$species,$help,$file);
$port = 3306;
$host='ens-genomics1';
$user='ensadmin';

GetOptions 
  (
   "dnadb_host=s" => \$dnadb_host,
   "dnadb_port=s" => \$dnadb_port,
   "dnadb_user=s" => \$dnadb_user,
   'dnadb_name=s' => \$dnadb_name,
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

#Core db not really important here but may pass core db as parameters to avoid errors?
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
     -dnadb_name => $dnadb_name,
     );
$efg_dba->dbc->db_handle;

print ":: Loading FeatureType Gene Links from:\t$file\n";

my $fta	= $efg_dba->get_FeatureTypeAdaptor;
my $dbentry_adaptor = $efg_dba->get_DBEntryAdaptor;

my $analysis_adaptor = $efg_dba->get_AnalysisAdaptor;
my $analysis = $analysis_adaptor->fetch_by_logic_name("Manual Antibody-Gene Annotation");

if(! defined $analysis){ throw "Analysis 'Manual Antibody-Gene Annotation' not found"; }

my $helper = Bio::EnsEMBL::Funcgen::Utils::Helper->new(no_log => 1); #default to STDOUT

my %releases;
open(FILE,$file);
<FILE>; #Throw away header

while(<FILE>){
  chomp;
  my ($feature_type, $gene_stable_id, $release) = split(/\t/);

  eval{	
    
    if(!defined($releases{$release})){
      #Try to create it...
      add_external_db($efg_dba,$species.'_core_Gene',$release,'EnsemblGene');
      $releases{$release}=1;
    }

    my $ft = $fta->fetch_by_name($feature_type);
    if(!$ft){ warn $feature_type." does not exist"; next; }
    
    my $gene_name = $helper->get_core_display_name_by_stable_id($efg_dba->dnadb(),$gene_stable_id,'Gene');
    
    my $dbentry = Bio::EnsEMBL::DBEntry->new(
					     -dbname             => $species.'_core_Gene',
					     -release            => $release,
					     -status             => 'KNOWNXREF',
					     #-display_label_linkable => 1,
					     -db_display_name    => 'EnsemblGene',
					     -type               => 'MISC', 
					     -primary_id         => $gene_stable_id,
					     -display_id         => $gene_name,
					     -info_type          => 'MISC',
					     -info_text          => 'GENE',
					     -linkage_annotation => $analysis->logic_name,
					     -analysis           => $analysis,
					     #-description        => 'Ensembl Gene associated to Feature Type',
					    );
    
    #1 as the last argument to ignore Gene release (edb version)
    $dbentry_adaptor->store( $dbentry, $ft->dbID, 'FeatureType');
  };

  if($@){ 
    warn $@;
  }
  else{
    #They may already have been present
    print "Succesfully stored Gene:FeatureType xref: $gene_stable_id:$feature_type\n";
  }
}

close FILE;
