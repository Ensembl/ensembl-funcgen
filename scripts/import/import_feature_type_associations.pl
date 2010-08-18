#!/usr/local/ensembl/bin/perl

=head1 NAME

ensembl-efg import_feature_type_associations.pl
  
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

my %fts;
open(FILE,"types/FeatureType_associations.txt");
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

