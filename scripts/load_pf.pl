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

=cut

use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::FeatureSet;
use Bio::EnsEMBL::Utils::Exception qw( throw );

$| =1;

my ($file, $ftype_name, $ctype_name, $pass, $line, $fset_id, %slice_cache);
my ($fset, $dbhost, $dbname, $cdbname, $help, $man, @p_features);
my $port = 3306;
my $anal_name = 'Nessie';


GetOptions (
	    "feature_type=s"   => \$ftype_name,
	    "file|f=s"         => \$file,
	    "cell_type=s"      => \$ctype_name,
	    "feature_set_id=i" => \$fset_id,
	    "pass=s"           => \$pass,
	    "port=s"           => \$port,
            "dbname=s"         => \$dbname,
	    "dbhost=s"         => \$dbhost,
            "cdbname=s"        => \$cdbname,
	    "analysis_name=s"  => \$anal_name,
	    "help|?"           => \$help,
	    "man|m"            => \$man,
	   );

pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;


### Set up adaptors and FeatureSet 

if(! $cdbname  || ! $dbname ){
  throw("You must provide a funcgen(-dbname) and a core(-cdbname) dbname");
}
throw("Must define your funcgen dbhost -dbhost") if ! $dbhost;
throw("Must supply an input file with -file") if ! $file;
throw("Must supply a password for your fungen db") if ! $pass;


my $cdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
					      -host => "ens-livemirror",
					      -dbname => $cdbname,
					      #-species => "homo_sapiens",
					      -user => "ensro",
					      -pass => "",
					      #	-group => 'funcgen',
					      -port => '3306',
					     );

my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
						      -host => $dbhost,
						      -dbname => $dbname,
						      #-species => "homo_sapiens",
						      -user => "ensadmin",
						      -pass => $pass,
						      -dnadb => $cdb,
						      -port => '3306',
						     );

#should check db's here


my $pfa = $db->get_PredictedFeatureAdaptor();
my $fset_adaptor = $db->get_FeatureSetAdaptor();

if($fset_id){
  $fset = $fset_adaptor->fetch_by_dbID($fset_id);


  throw("Could not retrieve FeatureSet with dbID $fset_id") if ! $fset;

  warn("You are loading PredictedFeature using a previously stored FeatureSet:\n".
       "\tCellType:\t".$fset->cell_type->name()."\n".
       "\tFeatureType:\t".$fset->feature_type->name()."\n");

  #should also check types and anal if the have been set
  #should add more on which experiment/s this is associated with and feature_set/dat_set_name when we have implemented it.
  #also ask continute question or use force_import flag
  #we could also do a count on the PFs in the set to make sure we're know we're adding to a populated set.
  #hard to repair if we do load ontop of another feature set


}elsif(! ($ftype_name && $ctype_name)){
  throw("Must provide a FeatureType and a CellType name to load your PredictedFeatures");
}else{
  my $anal =  $db->get_AnalysisAdaptor->fetch_by_logic_name($anal_name);
  my $ftype = $db->get_FeatureTypeAdaptor->fetch_by_name($ftype_name);
  my $ctype = $db->get_CellTypeAdaptor->fetch_by_name($ctype_name);


  throw("No valid CellType available for $ctype_name") if ! $ctype; 
  throw("No valid FeatureType available for $ftype_name") if ! $ftype;
  throw("No valid Analysis available for $anal_name") if ! $anal;
  
  $fset = Bio::EnsEMBL::Funcgen::FeatureSet->new
    (
     -CELL_TYPE => $ctype,
     -FEATURE_TYPE => $ftype,
     -ANALYSIS => $anal,
    );

  ($fset) = @{$fset_adaptor->store($fset)};

  warn("Generated FeatureSet\n");
}


open (FILE, $file) || die "Unable to open input\t$file";

warn("Loading PredictedFeatures from $file\n");

while ($line = <FILE>){

  chomp $line;
  my @tmp = split /\s+/, $line;
  my $start = $tmp[1];
  my $end = $tmp[2];
  my $score = $tmp[4];
  my $text = "enriched_site";
  my $chr = $tmp[0];

  $chr =~ s/chr//;

  #	print STDERR "$cdb->get_SliceAdaptor()->fetch_by_region(\'chromosome\', $chr);\n";
  

  if (! exists  $slice_cache{$chr}){
    $slice_cache{$chr} = $cdb->get_SliceAdaptor()->fetch_by_region('chromosome', $chr);
    throw("Could not generate slice for chromosome $chr") if ! $slice_cache{$chr};
  }

  my $pfeature = Bio::EnsEMBL::Funcgen::PredictedFeature->new
    (
     -SLICE         => $slice_cache{$chr},
     -START         => $start,
     -END           => $end,
     -STRAND        => 1,
     -DISPLAY_LABEL => $text,
     -SCORE         => $score,
     -FEATURE_SET   => $fset,
    );
  
  push @p_features, $pfeature;

}
my @chrs = keys(%slice_cache);
$pfa->store(@p_features);

warn("Loaded ".($.)." PredictedFeatures onto chromosomes @chrs\n");



__END__
