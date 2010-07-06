#!/usr/local/ensembl/bin/perl

=head1 NAME

ensembl-efg import_feature_type_gene_links.pl
  
=head1 SYNOPSIS

import_feature_type_gene_links.pl [options]

Options:

Mandatory

  -species|s       Species name

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
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );

my ($pass,$port,$host,$user,$dbname,$species,$help);
$port = 3306;
$host='ens-genomics1';
$user='ensadmin';

GetOptions (
	    "dbpass|p=s"   => \$pass,
	    "dbport|l=s"   => \$port,
	    "dbhost|h=s"   => \$host,
	    "dbuser|u=s"   => \$user,
	    "dbname=s"     => \$dbname,
	    "species|s=s"  => \$species, 
	    "help|?"       => \$help,
	   );
pod2usage(1) if $help;
pod2usage(0) if !defined($species);
pod2usage(0) if !defined($dbname);
pod2usage(0) if !defined($pass);

my $efg_dba = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
    (
     -host    => $host,
     -user    => $user,
     -dbname  => $dbname,
     -species => $species,  
     -port    => $port,
     -pass	  => $pass,
     );
my $fta	= $efg_dba->get_FeatureTypeAdaptor();     
my $dbentry_adaptor = $efg_dba->get_DBEntryAdaptor();

my $helper = Bio::EnsEMBL::Funcgen::Utils::Helper->new(
													   no_log => 1,#default to STDOUT
													  );

open(FILE,"types/".$species.".FeatureType_Genes.txt");
<FILE>; #Title
while(<FILE>){
	chomp;
	my ($feature_type, $gene_stable_id, $release) = split(/\t/);
		
	eval{	
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
			-linkage_annotation => 'ENSEMBL Manual Curation',
			-description        => 'ENSEMBL Gene associated to Feature Type',
		);

		#1 as the last argument to ignore Gene release
		$dbentry_adaptor->store( $dbentry, $ft->dbID, 'FeatureType');
	};
	if($@){ warn $@; };
}

close FILE;
