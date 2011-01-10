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

ensembl-efg convert_htilist_to_features.pl
  
=head1 SYNOPSIS

convert_hitlist_to_features.pl [options]

Options:

Mandatory


Optional


=head1 OPTIONS

=over 8

=item B<-name|n>

Mandatory:  Instance name for the data set, this is the directory where the native data files are located

=item B<-format|f>

Mandatory:  The format of the data files e.g. nimblegen

=over 8

=item B<-group|g>

Mandatory:  The name of the experimental group

=over 8

=item B<-data_root>

The root data dir containing native data and pipeline data, default = $ENV{'EFG_DATA'}

=over 8

=item B<-fasta>

Flag to turn on dumping of all probe_features in fasta format for the remapping pipeline

=item B<-norm>

Normalisation method, deafult is the Bioconductor vsn package which performs generalised log ratio transformations

=item B<-species|s>

Species name for the array.

=item B<-debug>

Turns on and defines the verbosity of debugging output, 1-3, default = 0 = off

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> takes a input redundant probe name fasta file and generates an NR probe dbID fasta file.

=cut


#add @INC stuff here, or leave to .bashrc/.efg?

BEGIN{
  if (! defined $ENV{'EFG_DATA'}) {
	if (-f "~/src/ensembl-functgenomics/scripts/.efg") {
	  system (". ~/src/ensembl-functgenomics/scripts/.efg");
	} else {
	  die ("This script requires the .efg file available from ensembl-functgenomics\n".
		   "Please source it before running this script\n");
	}
  }
}
	

#use Bio::EnsEMBL::Root; #Only used for rearrange see pdocs
#Roll own Root object to handle debug levels, logging, dumps etc.

### MODULES ###
use Getopt::Long;
#use Carp;#For dev only? cluck not exported by default Remove this and implement in Helper
use Pod::Usage;
#POSIX? File stuff
use File::Path;
#use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw (open_file run_system_cmd backup_file);
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::FeatureType;
use Bio::EnsEMBL::Funcgen::FeatureSet;
use Bio::EnsEMBL::Funcgen::DataSet;
use Bio::EnsEMBL::Funcgen::AnnotatedFeature;
use Bio::EnsEMBL::Analysis;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use strict;

$| = 1;							#autoflush
my ($pass, $dbname, $help, $man, $species);
my ($set_name, $data_version, $clobber, $rset);
my $user = "ensadmin";
my $host = 'localhost';
my $port = '3306';

#Definitely need some sort of Defs modules for each array?

$main::_debug_level = 0;
$main::_tee = 0;


#Use some sort of DBDefs for now, but need  to integrate with Register, and have put SQL into (E)FGAdaptor?
#Use ArrayDefs.pm module for some of these, class, vendor, format?
#ArrayDefs would also contain paths to data and vendor specific parse methods?

GetOptions (
			"pass|p=s"     => \$pass,
			"port=s"     => \$port,
			"host|h=s"     => \$host,
			"user|u=s"     => \$user,
			"dbname|d=s"   => \$dbname,
			"species=s"    => \$species,
			"help|?"       => \$help,
			"man|m"        => \$man,
			"set_name|n=s"   => \$set_name,
			"data_version|s=s" => \$data_version,
			'clobber'          => \$clobber,
		   );




pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;


my $cdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
											  -host => 'ensembldb.ensembl.org',
											  -user => 'anonymous',
											  -dbname => $species."_core_".$data_version,
											  -species => $species,
											 );



my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
													   -dbname => $dbname,
													   -port   => $port,
													   -pass   => $pass,
													   -host   => $host,
													   -user   => $user,
													   -dnadb  => $cdb,
													 );


my $fset_a = $db->get_FeatureSetAdaptor();
my $dset_a = $db->get_DataSetAdaptor();
my $rset_a = $db->get_ResultSetAdaptor();
my $anal_a = $db->get_AnalysisAdaptor();
my $slice_a = $db->get_SliceAdaptor();
my $afa = $db->get_AnnotatedFeatureAdaptor();
my $anal = $anal_a->fetch_by_logic_name('Parzen');

my @rsets = @{$rset_a->fetch_all_by_name_Analysis($set_name.'_IMPORT', $anal)};
if(scalar(@rsets) > 1){
  throw($set_name.' has more than one ResultSet with analysis '.$anal->logic_name());
}
else{
 $rset = $rsets[0];
}

throw("No ResultSet with name $set_name and analysis\t".$anal->logic_name()) if ! $rset;


my $fset = $fset_a->fetch_by_name($set_name);


if($fset && ! $clobber){
  throw("Found pre-existing FeatureSet:\t$set_name\nUse -clobber to overwrite");
}elsif($fset && $clobber){
  my $cs_id = $db->get_FGCoordSystemAdaptor->fetch_by_name('chromosome')->dbID();
  my $sql = 'DELETE from annotated_feature where feature_set_id='.$fset->dbID().' and coord_system_id='.$cs_id;
  #$db->dbc->do($sql) || throw('Failed to roll back annotated_features for feature_set_id'.$fset->dbID());
}else{#no fset
  $fset = Bio::EnsEMBL::Funcgen::FeatureSet->new(
												 -name         => $set_name,
												 -analysis     => $anal,
												 -feature_type => $rset->feature_type(),
												 -cell_type    => $rset->cell_type(),
												 -feature_class=> 'annotated',
												);

  ($fset) = @{$fset_a->store($fset)};
}

#my @tmp = @{$dset_a->fetch_all_by_FeatureSet($fset)};

if(! @{$dset_a->fetch_all_by_FeatureSet($fset)}){
  my $dset = Bio::EnsEMBL::Funcgen::DataSet->new(
												 -feature_set => $fset,
												 -result_set  => $rset,
												 -name        => $set_name,
												);

  $dset_a->store($dset);

}


#foreach my $slice(@{$slice_a->fetch_all('top_level')}){


#my @chrs = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y');
my @chrs = ('MT');

foreach my $chr (@chrs){

  my $cnt = 0;
  my $slice = $slice_a->fetch_by_region('chromosome', $chr);
  my @tmp = @{$rset->get_ResultFeatures_by_Slice($slice)};

  foreach my $rf(@tmp){

   	my $af = Bio::EnsEMBL::Funcgen::AnnotatedFeature->new(
														  -slice => $slice,
														  -start => $rf->start(),
														  -end   => $rf->end(),
														  -score => $rf->score(),
														  -feature_set => $fset,
														  -strand => 0,
														 );

	$afa->store($af);
	$cnt ++;
  }
  print "Loaded $cnt AnnotatedFeatures for $set_name on chromosome $chr\n";

}

1;
