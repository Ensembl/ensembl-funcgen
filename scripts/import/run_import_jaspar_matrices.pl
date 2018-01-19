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

ensembl-efg run_import_jaspar_matrices.pl
  
=head1 SYNOPSIS

run_import_jaspar_matrices.pl [options]

Options:

Mandatory

  -species|s       Species name

  -matrix_data_dir|d     Folder where Jaspar data is (both matrix_only.txt and matrix_list.txt)

  -dbpass|p        The password for the EFG DB

  -dbname          Defines the eFG dbname
  
Optional

  -file|f    File manually linking matrix to genes

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

Folder where the Jaspar matrix data is (matrix_only.txt and matrix_list.txt)

=item B<-pass|p>

Database password

=item B<-dbname>

Database name

=item B<-file|f>

File with links between matrix and genes
This file is text, each line comprising of: 
Jaspar_ID\tENSEMBL_GENE_1;ENSEMBL_GENE_2;..;ENSEMBL_GENE_N

By default, link between matrix and feature type is defined using internal core links
By providing this file internal links will be ignored.

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

=item B<-help>

Print a brief help message and exits.

=back

=head1 DESCRIPTION

B<This script> Imports matrices from Jaspar database files

By default infers the associated FeatureType from the protein associated to each matrix 
Feature Types must be loaded in the database and properly associated to genes.
Use the following scripts:
import_type.pl 
import_feature_type_associations.pl
import_feature_type_gene_links.pl
(or use CreateDB from the efg environment - see introduction_to_efg doc)

By default, a matrix is associated to feature types via internal core links using 
proteins associated to the matrix. A user can also provide a file manually linking
matrices to the genes of the species for the database where you're loading matrices

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

my ($dnadb_host,$dnadb_port,$dnadb_user);
my ($pass,$port,$host,$user,$dbname,$species,$matrix_data_dir,$file,$help);
$port = 3306;
$host='ens-genomics1';
$user='ensadmin';

GetOptions (
	    "dnadb_host=s" => \$dnadb_host,
	    "dnadb_port=s" => \$dnadb_port,
	    "dnadb_user=s" => \$dnadb_user,
	    "pass|p=s"    => \$pass,
	    "port|l=s"    => \$port,
	    "host|h=s"    => \$host,
	    "user|u=s"    => \$user,
	    "file|f=s"    => \$file,
	    "dbname=s"      => \$dbname,
	    "species|s=s"   => \$species, 
	    "matrix_data_dir=s"    => \$matrix_data_dir,
	    "help|?"        => \$help,
	   );
pod2usage(1) if $help;
pod2usage(0) if !defined($species);
pod2usage(0) if !defined($matrix_data_dir);
pod2usage(0) if !defined($dbname);
#pod2usage(0) if !defined($pass);

my $efg_dba = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
    (
     -host    => $host,
     -user    => $user,
     -dbname  => $dbname,
     -dnadb_host => $dnadb_host,
     -dnadb_port => $dnadb_port,
     -dnadb_user => $dnadb_user,
     -species => $species,  
     -port    => $port,
     -pass    => $pass,
     );

if(!(-e $matrix_data_dir."/matrix_only.txt")){
  throw "Expecting ".$matrix_data_dir."/matrix_only.txt"; 
}

if(!(-e $matrix_data_dir."/matrix_list.txt")){
  throw "Expecting ".$matrix_data_dir."/matrix_list.txt"; 
}

#Cache matrix data
my %matrix_data;
open(FILE,$matrix_data_dir."/matrix_only.txt");
while(<FILE>){
  if(/^>(\S+)\s*/){
    my $matrix_id = $1;
    my $matrix = <FILE>.<FILE>.<FILE>.<FILE>;
    $matrix_data{$matrix_id} = $matrix;
  }
}
close FILE;


my $aa	= $efg_dba->get_AnalysisAdaptor();
my $analysis = $aa->fetch_by_logic_name('Jaspar');
if(!$analysis){ throw "An Analysis named Jaspar needs to be in the EFG database"; }

my $ga = $efg_dba->dnadb->get_GeneAdaptor();
my $fta = $efg_dba->get_FeatureTypeAdaptor();
my $bma	= $efg_dba->get_BindingMatrixAdaptor();
my $dbea = $efg_dba->get_DBEntryAdaptor();

#We will only upload data for TFs so only search for those.
my @tfs = @{$fta->fetch_all_by_class('Transcription Factor')};
push @tfs, @{$fta->fetch_all_by_class('Insulator')};
push @tfs, @{$fta->fetch_all_by_class('Transcription Factor Complex')};

#Just caching data... assume ft->name is unique at least among TFs
my %feature_type;
map { $feature_type{$_->name} = $_; } @tfs;

my %feature_type_genes;
#Get the genes associated to the TF and associated TFs (relevant for complexes)
foreach my $tf (@tfs){
  my @assoc_tfs =  @{$fta->fetch_all_by_association($tf)};
  push @assoc_tfs, $tf;
  my @genes;
  foreach my $assoc_tf (@assoc_tfs){
    #Assume all dbentries associated to feature type are genes?? otherwise force checking
    my @db_entries = @{$dbea->fetch_all_by_FeatureType($assoc_tf)};
    next if(scalar(@db_entries)==0);
    if(scalar(@db_entries)>1){ 
      warn "Dubious association to gene in TF ".$tf->name." in ".$_->name; 
      next; 
    }
    push @genes, $db_entries[0]->primary_id;
  }
  next if(scalar(@genes)==0);
  $feature_type_genes{$tf->name} = \@genes;
}


my %matrix_genes;
if($file){
  #Get the links from the file
  open(FILE,$file);
  while(<FILE>){
    chomp;
    my ($matrix,$gene_list) = split(/\t/);
    my @gene_ids = split(/;/,$gene_list);
    $matrix_genes{$matrix} = \@gene_ids;
  }
  close FILE;

} else {
  #Run through the matrix list and get core links  
  open(FILE,$matrix_data_dir."/matrix_list.txt");
  while(<FILE>){
    chomp;
    my ($matrix,$ic,$name,$family,$proplist) = split(/\t/);
    $proplist =~ /acc\s+\"(\S+)\"\s+/;
    next if(!$1 || ($1 eq '-'));
    my @acclist = split(/,/,$1);
    my @gene_ids;
    foreach my $acc (@acclist){
      my @genes = @{$efg_dba->dnadb->get_GeneAdaptor()->fetch_all_by_external_name($acc)};
      if(scalar(@genes)>1){ warn $acc." associated to more than one gene... ambiguous reference ignored"; }
      if(scalar(@genes)==1){ push @gene_ids, $genes[0]->stable_id; }
    }
    next if(scalar(@gene_ids)==0);
    $matrix_genes{$matrix} = \@gene_ids;
  }
  close FILE;
  
}

foreach my $matrix (sort keys %matrix_genes){
  my @m_genes = @{$matrix_genes{$matrix}};
  foreach my $ft (sort keys %feature_type_genes){
    my @ft_genes = @{$feature_type_genes{$ft}};
    #Only import if ALL genes are the same between matrix and TF
    next if(scalar(@m_genes)!=scalar(@ft_genes));
    #I have the feeling there is an easier way of doing this!!
    my $not_found=0; 
    foreach my $g1 (@m_genes){ 
      my $found=0;
      foreach my $g2 (@ft_genes){ 
	if($g1 eq $g2){ $found=1; last; } 
      }
      if($found==0){ $not_found=1; last; }
    }
    next if($not_found==1);

    #Check if it is not already imported...
    if(@{$bma->fetch_all_by_name_FeatureType($matrix,$feature_type{$ft})}){ 
      warn "Matrix $matrix already exists in the database associated to $ft. Skipping import"; 
    } else {
      #If not, then import
      my $bm = Bio::EnsEMBL::Funcgen::BindingMatrix->new(
							 -name         => $matrix,
							 -analysis     => $analysis,
							 -feature_type => $feature_type{$ft},
							 -description  => 'Jaspar Matrix',
							 -frequencies  => $matrix_data{$matrix}
							);
      
      ($bm) = $bma->store($bm);
      if(!$bm){ 
	warn "Failed import of $matrix, associated to $ft"; 
      }
    }
  }
}
