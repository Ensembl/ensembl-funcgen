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

BEGIN{
	if(! defined $ENV{'EFG_DATA'}){
		if(-f "$ENV{'HOME'}/src/ensembl-functgenomics/scripts/.efg"){
			system (". ~/src/ensembl-functgenomics/scripts/.efg");
		}else{
			die ("This script requires ensembl-functgenomics/scripts/.efg\n".
				 "Please source it before running this script\n");
		}
	}
}

use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Funcgen::Helper;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;#

use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(species_chr_num open_file);
use Bio::EnsEMBL::Utils::ConfigRegistry;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::AnnotatedFeature;
use Bio::EnsEMBL::Funcgen::RegulatoryFeature;
use Bio::EnsEMBL::Funcgen::DataSet;
use File::Path qw(mkpath);
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(mean);
use Devel::Size qw(size total_size);
use Bio::EnsEMBL::Funcgen::ResultFeature;
use Tie::File;
#use Bio::MAGE::XMLUtils;
use Bio::MAGE::Experiment::Experiment;

$| = 1;
use strict;
my $reg = "Bio::EnsEMBL::Registry";


#my @rfs;
#for my $i(0..10000000){
#  push @rfs, Bio::EnsEMBL::Funcgen::ResultFeature->new_fast([0, 50, 2]);

  #acess methods but don't print 
#  $rfs[$i]->score();
#    $rfs[$i]->start();
#  $rfs[$i]->end();
  
#}

#print 'Size is '.size(\@rfs)."\n";
#print 'Total Size is '.total_size(\@rfs)."\n";

#exit;


my $cdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
  (
   -host => 'ens-staging',#livemirror',#"ensdb-1-12",
   #-host => "ensdb-archive",
   -dbname => 'homo_sapiens_core_46_36h',
   #-dbname => 'homo_sapiens_funcgen_45_36g',
   -species => "Homo_sapiens",
   -user => "ensro",
   #-user => 'anonymous',
  #-pass => "",
   #-group => 'core',
   #-port => '3304',
  );

my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
  (
   -host => 'ens-genomics1',#livemirror',#"ensdb-1-12",
   #-host => "ensdb-archive",
   #-dbname => 'homo_sapiens_funcgen_46_36h',
   -dbname => 'homo_sapiens_funcgen_47_36i',
   -species => "Homo_sapiens",
   -dnadb => $cdb,
   -user => "ensadmin",
   #-user => 'anonymous',
   -pass => "ensembl",
   #-group => 'core',
  # -port => '3307',
  );





#my $afa = $db->get_AnnotatedFeatureAdaptor();
#my $pfa = $db->get_PredictedFeatureAdaptor();
my $slice_a = $db->get_SliceAdaptor();
my $slice = $slice_a->fetch_by_region('chromosome', 14);#, 23083751, 23084751);#23118316);





my $fset = $db->get_FeatureSetAdaptor->fetch_by_name('RegulatoryFeatures_v45');

my (%ftype_hash, $stable_id, $sql);

#reg feature hack
my %reg_class_regexs = (
						#'1....(10|01).'  => 'Gene end associated', 
						#'1...1...'        => 'Promoter associated',#orig
						'1...1.....'        => 'Promoter associated',
						'1.0.001...' => 'Non-gene associated',
						'11..01....' => 'Gene associated',
					   );



#omit TSS and TES from here?
my @reg_feature_attrs = ('DNase1', 'CTCF', 'H4K20me3', 'H3K27me3', 
						 'H3K36me3', 'H3K4me3', 'H3K79me3', 'H3K9me3', 'TSS Proximal', 'TES Proximal'); 

#have to do it one by one other wise we'll run out of memory

my $ft_adaptor = $db->get_FeatureTypeAdaptor();

my $feat_adaptor = $fset->get_FeatureAdaptor();

foreach my $dbID(@{$feat_adaptor->list_dbIDs()}){


  my $regf = $feat_adaptor->fetch_by_dbID($dbID);

warn " got regf $regf with dbID $dbID";

  #this is defined tho as we're dynamically setting the type
  if(! defined $regf->feature_type() && $regf->display_label =~/[01]/){
	my $ftype;

   #my @vector = split//, $display_label;
		  
	  #foreach my $i(0..7){#$#vector){
	#	push @$reg_attrs, $reg_feature_attrs[$i] if $vector[$i];
	#	  }
		
		  
	foreach my $regex(keys %reg_class_regexs){
	  
	  if($regf->display_label() =~ /$regex/){
		
		#warn "$vector matches ".$reg_class_regexs{$regex}."\t$regex\n";
		
		throw('Found non-mutually exclusive regexs') if $ftype;
		$ftype = $reg_class_regexs{$regex};
	  }
	  
	}
	
	#  undef $display_label;
	$ftype ||= 'Unclassified';
	$ftype_hash{$ftype} = $ft_adaptor->fetch_by_name($ftype) if (! exists $ftype_hash{$ftype});
	$ftype = $ftype_hash{$ftype};

	($stable_id = $regf->stable_id()) =~ s/ENSR0+//;

	$sql = 'UPDATE regulatory_feature set feature_type_id='.$ftype->dbID().', stable_id="'.$stable_id.'" where regulatory_feature_id='.$regf->dbID();

	$db->dbc->do($sql);
  }

}


