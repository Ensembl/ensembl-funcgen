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

remap_features.pl

=head1 SYNOPSIS

this script will ...

=head1 DESCRIPTION

=cut

use strict;
use warnings;

use Data::Dumper;
use Getopt::Long;

use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw (open_file);
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Funcgen::PredictedFeature;


my ($pass, $dbname, $odata_version, $ndata_version, $man, $help, $species, $clobber);
my ($old_assembly, $new_assembly, @chrs, $set_name, $dump_features, $write_features);

my $out_dir = ".";
my $host = 'ens-genomics2';
my $user = 'ensadmin';
my $port = 3306;

GetOptions (
			"pass|p=s"         => \$pass,
			"port=s"           => \$port,
			"host|h=s"         => \$host,
			"user|u=s"         => \$user,
			"dbname|d=s"       => \$dbname,
			"species=s"        => \$species,
			"help|?"           => \$help,
			"man|m"            => \$man,
			"set_name=s"     => \$set_name,
			"old_data_version=s" => \$odata_version,#of core db which contains old to new assembly mapping
			"new_data_version=s" => \$ndata_version,
			"old_assembly|o=s" => \$old_assembly,
			"new_assembly|n=s" => \$new_assembly,
			'clobber'          => \$clobber,
			"write_features"   => \$write_features,
			"dump_features"    => \$dump_features,
			"out_dir=s"        => \$out_dir,
			"chr_name=s"       => \@chrs,
		   );


if(! @chrs){
  #need to add alis for encode regions here
  print "No chromosome specified, running wiht default all chromosomes\n";
  @chrs = (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y');
}else{
  @chrs = split/,/, join(',', @chrs);
}

#set default here and warn, rather than above
#check all before throwing

throw('Must provide a dbname') if ! defined $dbname;
throw('Must provide a password') if ! defined $pass;
throw('Must provide a species') if ! defined $species;
throw('Must provide a old_data_version') if ! defined $odata_version;
throw('Must provide a new_data_version') if ! defined $ndata_version;
throw('Must provide an old_assembly to map from') if ! defined $old_assembly;
throw('Must provide an new_assembly to remap to') if ! defined $new_assembly;
throw('Must provide a feature_set name') if ! defined $set_name;
if (! ($write_features || $dump_features)) {
    print "No output type specified turning on dump_features\n";
    $dump_features = 1;
}


$| = 1;


# connect to databases and get adaptors
use Bio::EnsEMBL::Registry;

#can't load from registry as this will get the latest dnadb
#we may want to map from an earlier DB
my $ocdbname = $species.'_core_'.$odata_version;
my $dnadb = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
												-host   => 'ensembldb.ensembl.org',
												-port   => 3306,
												-user   => 'anonymous',
												#-pass   => $dbpass,
												-dbname => $ocdbname,
											   );




my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
  (
   -host   => $host,
   -port   => $port,
   -user   => $user,
   -pass   => $pass,
   -dbname => $dbname,
   -dnadb  => $dnadb
  );

my $source_sa = $db->get_SliceAdaptor();
my $fsa = $db->get_FeatureSetAdaptor();
my $source_pfa = $db->get_PredictedFeatureAdaptor();

my $ncdbname = $species.'_core_'.$ndata_version;
$dnadb =  Bio::EnsEMBL::DBSQL::DBAdaptor->new(
												-host   => 'ensembldb.ensembl.org',
												-port   => 3306,
												-user   => 'anonymous',
												#-pass   => $dbpass,
												-dbname => $ncdbname,
											   );

#reset the dnadb to the one containing the assembly mapping
$db->dnadb($dnadb);
my $mapper_sa =  $db->get_SliceAdaptor();
my $mapper_pfa = $db->get_PredictedFeatureAdaptor();
my $aa = $db->get_AnalysisAdaptor();
my $dsa = $db->get_DataSetAdaptor();
my $out_file = $out_dir."/${set_name}.${new_assembly}";

#should do some extensive testing of CS's here to make sure mapping is available
#get meta container and check assembly.mapping
#check the default assembly in the mapper db is new_assembly

my $fset = $fsa->fetch_by_name($set_name);
throw("FeatureSet name $set_name does not exist in $dbname") if ! defined $fset;

my $fh = open_file($out_file, '>');


if($clobber && $write_features){
  my $fg_cs = $db->get_FGCoordSystemAdaptor->fetch_by_name('chromosome', $new_assembly);
  my $sql = 'delete from predicted_feature where feature_set_id='.$fset->dbID().' and coord_system_id='.$fg_cs->dbID();
  
  print "Deleting PredictedFeatures from $set_name on $new_assembly\n";
  $db->dbc->do($sql) or throw("Failed to delete PredictedFeatures from $set_name on $new_assembly");
}



foreach my $chr (@chrs){
  my (@pfs);
  my $cnt = 0;
  print "Projecting features on chromosome $chr\n";
  
  my $slice = $source_sa->fetch_by_region('chromosome', $chr, undef, undef,
										  undef, $old_assembly);
  throw("Assembly version in database does not match given option!")
    if ($slice->coord_system()->version() ne $old_assembly);



  #This is a slice based on the seq_regions_ids in the new DB, 
  #not the seq_region_ids from the original import DB which can/will be different



  #we to query the DB using the original import DB
  #then reset the dnadb to the db with the mapping
  #the for each feature recreate it using the new dnadb slice with the old assembly
  #then project to new assemlby
 
  my $mapper_slice_oldasm = $mapper_sa->fetch_by_region('chromosome', $chr, undef, undef,
												 undef, $old_assembly);
  my $mapper_slice_newasm = $mapper_sa->fetch_by_region('chromosome', $chr, undef, undef,
												 undef, $new_assembly);
  


  print "Fetching PredictedFeatures from:\t$ocdbname\n";
  my $fts = $source_pfa->fetch_all_by_Slice_FeatureSet($slice, $fset);
  
  print "Found ".scalar(@$fts)." PredictedFeatures on $old_assembly\n";

  foreach my $ft (@$fts) {
    
    #print $fh join("\t", $fset_name, $ft->seq_region_name, $ft->start, $ft->end, $ft->strand), "\n";
  

    # create a new feature on the old assembly using the mapper db
    my $pf = Bio::EnsEMBL::Funcgen::PredictedFeature->new
	  (
	   -slice         => $mapper_slice_oldasm,
	   -start         => $ft->start,
	   -end           => $ft->end,
	   -strand        => $ft->strand,
	   -display_label => $ft->display_label,
	   -score         => $ft->score,
	   -feature_set   => $fset,
	  );
    
    #project feature to new assembly
    my @segments = @{ $pf->feature_Slice->project('chromosome', $new_assembly) };
    #print Dumper @segments;
    
    # do some sanity checks on the projection results:
    # discard the projected feature if
    #   1. it doesn't project at all (no segments returned)
    #   2. the projection is fragmented (more than one segment)
    #   3. the projection doesn't have the same length as the original
    #      feature
 

	#add logging of unprojected features here
   
    # this tests for (1) and (2)
    if (scalar(@segments) == 0) {
        print "Feature doesn't project!\n";
        next;
    } elsif (scalar(@segments) > 1) {
        print "Projection of feature is fragmented!\n";
        next;
    }
    
    # test (3)
    my $proj = $segments[0]->to_Slice;
    if ($pf->length != $proj->length) {
	  warn("Projection and original feature differ in length!");
    }
    
    # everything looks fine, so adjust the coords of your feature
    $pf->start($proj->start);
    $pf->end($proj->end);
    $pf->slice($mapper_slice_newasm);
    
    print $fh join("\t", $pf->display_label, $pf->seq_region_name, $pf->start, $pf->end, $pf->strand), "\n" if $dump_features;
    $cnt ++;
    push @pfs, $pf if $write_features;
    
  }

  
  print "Projected $cnt PredictedFeatures for chromosome $chr to $new_assembly\n";

  if($write_features){
	print "Loading PredictedFeatures\n";




	$mapper_pfa->store(@pfs);
  }
}



#    my $analysis = $aa->fetch_by_logic_name('projected_feature');
#    if (! defined $analysis) {
#        
#        # create an analysis for the type of feature you wish to store
#        $analysis = new Bio::EnsEMBL::Analysis
#            (
#             -logic_name      => 'projected_feature',
#             -program_file    => 'remap_features.pl',
#             -parameters      => '-o NCBI36 -n NCBI35',
#             -description     => 'predicted feature remapped from NCBI36 to MCBI35.',
#             #-display_label   => 'enriched_site_v35',
#             -displayable     => 1
#             );
#        $analysis = $aa->store($analysis);
#    }
#
#    ### FeatureSet
#    my $new_fset_name = $fset_name."_".$opts{n};
#    my $new_fset = $fsa->fetch_all_by_name($new_fset_name);
#    if (! @{$new_fset}) {
#        # add new feature set
#        $new_fset = Bio::EnsEMBL::Funcgen::FeatureSet->new
#            (
#             -analysis => $analysis,
#             -feature_type => $fset->feature_type,
#             -cell_type => $fset->cell_type,
#             -name => $new_fset_name
#             );
#        $new_fset = $fsa->store($new_fset);
#    } else {
#        throw("More than one feature_set is currently not supported.")
#            if (scalar(@{$new_fset}) > 1);
#    }
#    $new_fset = $new_fset->[0];
# 
###
### Need to make sure that feature_set is linked to result_set
### 
#    ### DataSet
#    my $new_dset = $dsa->fetch_all_by_FeatureSet($new_fset);
#    my $add_dset = 0;
#    if (! @{$new_dset}) {
#        
#        warn("Need to add new data set(s)!");
#        foreach my $rset (@{$self->result_sets()}) {
#            # add new data set
#            $dset = Bio::EnsEMBL::Funcgen::DataSet->new
#                (
#                 -result_set => $rset,
#                 -feature_set => $new_fset,
#                 #-name => $fset_name
#                 );
#            $dsa->store($dset);
#        }
#
#    } else {
#        
#        
#    }
###
###
###
    
    
#    my $pf = $pfa->fetch_all_by_Slice_FeatureSet($slice_newasm, $new_fset);
#    if (@$pf) {
#        
#        warn("Not storiing features. Slice ".
#              join(':', $slice_newasm->seq_region_name,$slice_newasm->start,$slice_newasm->end).
#              " already contains ".scalar(@$pf)." predicted features of feature set ".
#              $new_fset->dbID.".");
#        
#    } else {
#
        # store the features
#        $pfa->store(@pf);
    
#    }
#}


1;
