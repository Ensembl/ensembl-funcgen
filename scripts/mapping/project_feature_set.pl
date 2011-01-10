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

project_feature_set.pl - Projects a feature set to a new genome assembly.

=head1 SYNOPSIS

project_feature_set.pl [arguments]

Required arguments:

  --dbname, db_name=NAME              database name NAME
  --host, --dbhost, --db_host=HOST    database host HOST
  --port, --dbport, --db_port=PORT    database port PORT
  --user, --dbuser, --db_user=USER    database username USER
  --pass, --dbpass, --db_pass=PASS    database passwort PASS

Optional arguments:

  --conffile, --conf=FILE             read parameters from FILE
                                      (default: conf/Conversion.ini)

  --logfile, --log=FILE               log to FILE (default: *STDOUT)
  --logpath=PATH                      write logfile to PATH (default: .)
  --logappend, --log_append           append to logfile (default: truncate)
  --loglevel=LEVEL                    define log level (default: INFO)

  --is_component, --is-component      script is called from a wrapper script

  -i, --interactive                   run script interactively (default: true)
  -n, --dry_run, --dry                don't write results to database
  -h, --help, -?                      print help (this message)

=head1 DESCRIPTION

This script projects features from an old assembly to a new assembly

=cut

use strict;
use warnings;
#no warnings 'uninitialized';

use FindBin qw($Bin);
use Bio::EnsEMBL::Utils::ConfParser;
use Bio::EnsEMBL::Utils::Logger;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw( throw );

use Bio::EnsEMBL::Funcgen::Utils::Helper;#replace logger or inherit from logger?

#This is quite useful, lists params defined
#Summarises num wanring runtime, complete time and memusage at end.
# parse configuration and commandline arguments
my $conf = new Bio::EnsEMBL::Utils::ConfParser(
  -SERVERROOT => "$Bin/../../..",
  -DEFAULT_CONF => ""
);

$conf->parse_options
  (
   'host=s' => 1,
   'port=n' => 0,
   'user=s' => 1,
   'pass=s' => 0,
   'dbname=s' => 1,
   'old_cdbname=s' => 1,
   #'old_cdbport=n' => 0,
   #'old_cdbuser=s' => 0,
   #'old_cdbpass=s' => 0,
   #'old_cdbhost=s' => 0,
   'cdbname=s' => 0,
   'cdbport=n' => 0,
   'cdbuser=s' => 0,
   'cdbpass=s' => 0,
   'cdbhost=s' => 0,
   'clobber' => 0,
   'slice=s' => 0,
   'species=s' => 0,
   'ignore_length' => 0,
   #'force_store' => 1,
   'feature_set=s' => 1,
   'old_assembly=s' => 1,
   'new_assembly=s' => 1,
   'coord_system=s' => 0,
   'associations'      => 0,
  );

#$main::_no_log = 1;
#$main::_tee    = 1;#Isn't this set by default if no log i set?
my $helper = new Bio::EnsEMBL::Funcgen::Utils::Helper;

#assign assemblies for regex embedding
my $new_assembly = $conf->param('new_assembly');
my $old_assembly = $conf->param('old_assembly');
my $cs_level     = $conf->param('coord_system') || 'chromosome';


# get log filehandle and print heading and parameters to logfile
my $logger = new Bio::EnsEMBL::Utils::Logger(
  -LOGFILE    => $conf->param('logfile'),
  -LOGPATH    => $conf->param('logpath'),
  -LOGAPPEND  => $conf->param('logappend'),
  -VERBOSE    => $conf->param('verbose'),
);

# initialise log
$logger->init_log($conf->list_param_values);

# connect to database and get adaptors

#We only need the old DB to get access to the old toplevel
#So can use ensembldb
#This may include contigs which have now been merged into a chromosome
#and so would not be easily accessable via the new DB
#we would have to somehow get all non toplevel regions which are not 
#included in the new assembly.
#But they would not be present in the new DB
#We don't actually have this cross level mapping anyway, so can remove
#And only use new DB!

#Need to expose old_db params here!

my $old_cdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
  (
   -host   => 'ensdb-archive',
   -user   => 'ensro',
   -dbname => $conf->param('old_cdbname'),
   -group  => 'core',
   -port   => 5304,
  );

my ($new_cdb);

if(defined $conf->param('cdbname') || defined $conf->param('cdbhost') || defined $conf->param('cdbport') ){

  #This is the one with the mapping path between the assemblies
 
  $new_cdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
	(
	 -host   => $conf->param('cdbhost') || $conf->param('host'),
	 -port   => $conf->param('cdbport') || $conf->param('port'),
	 -user   => $conf->param('cdbuser') || $conf->param('user'),
	 -pass   => $conf->param('cdbpass') || $conf->param('pass'),
	 -dbname => $conf->param('cdbname'),
	 -group  => 'core',
	);
}elsif(! defined  $conf->param('species')){

  ##????????????????????????
  #This is now handled by the DBAdaptor and is not necessary if we have the latin species name in meta
  throw('Must provide a species name if no core DB has been defined');
}

my $efg_db = new Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor
  (
   -host   => $conf->param('host'),
   -port   => $conf->param('port'),
   -user   => $conf->param('user'),
   -pass   => $conf->param('pass'),
   -dbname => $conf->param('dbname'),
   -species =>  $conf->param('species'),
   -group  => 'funcgen',
   -dnadb  => $new_cdb,
  );



#Test connections
$efg_db->dbc->db_handle;
$new_cdb = $efg_db->dnadb;#set new_cdb if using default
$efg_db->dnadb->dbc->db_handle;
$old_cdb->dbc->db_handle;


#This is required for clobber
#As we need to project the slices to 
#see whether we have any previous projected features
#How will this handle two sequences which have been merged?
#We may get one set of features overwriting another
#We would have to do clobber at start before we project anything
#Actually NT are not mapped if they are merged so not so useful
#but may change?

#THis is wrong and we don't even need this as we can use the project method on slices from the new DB.
#This will obviosu not work for sr's which have disappeared between releases, but they are
#not caught by mapping anyway.

#my $ap = Bio::EnsEMBL::Utils::AssemblyProjector->new
#  (
#   -OLD_ASSEMBLY    =>  $old_assembly,
#   -NEW_ASSEMBLY    =>  $new_assembly,
#   -CHECK_LENGTH    => 0,
#   -MERGE_FRAGMENTS => 1,
#   -ADAPTOR         => $new_cdb,
#  );




my $old_sa = $old_cdb->get_SliceAdaptor;
my $new_sa = $new_cdb->get_SliceAdaptor;
my $fset_adaptor = $efg_db->get_FeatureSetAdaptor();
my $fset = $fset_adaptor->fetch_by_name($conf->param('feature_set'));	
								   
if(! defined $fset){
  throw("Cannot findFeatureSet:\t".$conf->param('feature_set')."\n");
}

my $set_feat_adaptor = $fset->get_FeatureAdaptor;
my $feat_class       = ucfirst($fset->type).'Feature';


#Do this slice wise to avoid possibility of running out of memmory for large sets
#or we could pull everythign back based on the feature_set?
my (@slices);#, $old_assembled_cs);

if($conf->param('slice')){
  @slices = ($old_sa->fetch_by_name($conf->param('slice')));

  warn('Need to test vs old assembly here');
}
else{
  #we can get toplevel for old assembly so we have to do for all?
  #Just do for chromosome now.
  #warn 'Hardcoded for just chromosome slices on old assembly due to no toplevel';
  #This will break for assemblies which does not have chromosome as top assembled level
  #These maybe scaffold, gene_scaffold, supercontig etc...
  #Currently mapping is only done for regions which are present in both assemblies
  #So NTs which are incorporated into chrs are not handled.


  #We just need to check that the default version is correct
  my $old_cs = $old_cdb->get_coordSystemAdaptor->fetch_by_name($cs_level);
  
  if($old_cs->version ne $old_assembly){
	warn "we need to implement old_assembled_level";
	warn "Need to change assumption that old level is default level";
	die("Old assembly($old_assembly) does not match default version in ".$conf->param('old_cdbname').":\t".$old_cs->version);
  }

  @slices = @{$old_sa->fetch_all($cs_level, $old_assembly)};
  #@slices = @{$old_sa->fetch_all('toplevel', undef, 1)};#inc non-ref
  #We can't do this on toplevel as projection only works for assembled levels
  #So long as other levels are not versioned, the features should persist.

  



  #We actually need to get the default assembled level with the given versionb
  #Normally chromosome but can be scaffold!
  #This also has to come from the old DB as we may want to remap from scaffold to chromosome?
  #This is why we need the old DB, to get the old toplevel which is not maynot be accessible in the new DB?
  #Compenents must be present for the mapping, but no easy way of accessing them.
  #my $old_assembled_cs;
  #No fetch_by_version_method??
  #foreach my $cs($old_cdb->get_CoordSystemAdaptor->fetch_all){
  #	if($cs->version eq $old_assembly){
  #	  $old_assembled_cs = $cs;
  #	  last;
  #	}
  #  }
  #if(! $old_assembled_cs){
  #	die("Could not find an old CoordSystem with version:\t$old_assembly");
  #  }
  #@slices = @{$sa->fetch_all($assembl_cs->name, $old_assembly)};

  if(! @slices){
	die("There are no slices available for the old assembly $old_assembly. Maybe it is not present in the DB?");
  }
}
	
#should check if we can get coord_system on new assembly here, or will this be caught by mapper?


$logger->info("Projecting features on ".scalar(@slices)." slices from assembly $old_assembly to assembly $new_assembly\n");

my ($new_slice, $slice_name, $failed_cnt, $wrong_length_cnt, $old_slice, $multi_segment_cnt);
my ($new_fslice, $old_fslice, $no_projection_cnt, $stored_cnt, @old_feats, @new_feats);
my $total_failed        = 0;
my $total_wrong_length  = 0;
my $total_stored        = 0;
my $total_multi_segment = 0;
my $total_no_proj       = 0;
my $length_txt = ($conf->param('ignore_length')) ? 'features were stored despite mismatched lengths in the new assembly' : 'features were skipped due to mismatched lengths in the new assembly';



#This will not currently work as the seq_region_id are most likely going to be different between the current cdb
#and the db on which the features were loaded.
#unless we are already mapping for the correct seq_region_id for the assemmbly
#This is failing because the coord systems are not equal
#can we get around this by resetting the dnadb?
#this will most likely cause problems as the seq_region_ids won't be consistent
#we could simply use the old DB to retrieve them and then set the new DB when we are storing
#this should be fine so long as we are changing the slice?

#Do check/clobber first just incase 
#we have merged slices which may overwrite each other
#if we rollback and store on a slice by slice basis

foreach my $old_slice(@slices){
  
  #redefine using new DB
  my $old_name = $old_slice->name;

  my $old_slice = $new_sa->fetch_by_name($old_name);

  if(! $old_slice){
	warn "Old seq_region is not present in new DB:\t$old_name\n";
	next;
  }

  my @segments = @{$old_slice->project($cs_level, $new_assembly)};

  #Now for each new segment do some rollback

  foreach my $segment(@segments){

	my $new_slice = $segment->to_Slice;

	#This should automatically use the correct seq_region_cache
	@new_feats = @{$fset->get_Features_by_Slice($new_slice)};
	
	if(@new_feats){
	  
	  if($conf->param('clobber')){
		#we actually want to force here as this will try and protect the old features
		#but we are only rolling back the new features
		$helper->rollback_FeatureSet($fset, 1 , $new_slice);
	  }
	  else{
		die(scalar(@new_feats).' projected features already exist for '.$old_slice->name.' > '.$new_slice->name."\nMaybe you want to set -clobber?\n");
	  }
	}
  }
}


#Now project the features
print "\nProjecting features...\n";
my $total_old_feats = 0;
my $dbentry_adaptor = $efg_db->get_DBEntryAdaptor();


foreach my $slice(@slices){
   
  @new_feats = ();

  #All we need to do now is regenrate the old slice in the newer dnadb which has the mapping
  #And then project them
  my $old_slice_name  = $slice->name;
  my $old_slice       = $new_sa->fetch_by_name($old_slice_name);
  (my $new_slice_name = $old_slice_name) =~ s/$old_assembly/$new_assembly/;
  my $new_slice = $new_sa->fetch_by_region($cs_level, $old_slice->seq_region_name);#No we can't do this here as we are assuming
  #That all projections will be to same slice

  @old_feats = @{$fset->get_Features_by_Slice($old_slice)};
  my $num_old_feats  = scalar(@old_feats);
  $total_old_feats  += $num_old_feats;
  $failed_cnt        = 0;
  $wrong_length_cnt  = 0;
  $stored_cnt        = 0;
  $multi_segment_cnt = 0;
  $no_projection_cnt = 0;


  #Project each feature

  if (! @old_feats){
	print "No features found for slice:\t".$old_slice->name."\n";
	next;
  }

  print "\nProjecting $num_old_feats features for slice:\t".$old_slice->name."\n";

  foreach my $feat(@old_feats){
		
	my @segments = @{$feat->project($cs_level, $new_assembly)};
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
	  #print "Feature doesn't project!\n";
	  $failed_cnt++;
	  $no_projection_cnt++;
	  next;
    } elsif (scalar(@segments) > 1) {
	  $failed_cnt++;
	  $multi_segment_cnt++;
	  next;
    }
    
    # test (3)
    my $proj_slice = $segments[0]->to_Slice;

    if ($feat->length != $proj_slice->length) {
	  $wrong_length_cnt++;

	  if(! $conf->param('ignore_length')){
		$failed_cnt++;
		next;
	  }
    }
    
    # everything looks fine, so adjust the coords of your feature
	#Have to generate new_slice here as we are not sure it is going to be 
	#on the same slice as the old assembly
	my $new_full_slice = $new_sa->fetch_by_region($cs_level, $proj_slice->seq_region_name);

    $feat->start($proj_slice->start);
    $feat->end($proj_slice->end);
    $feat->slice($new_full_slice);
    #do we need to deal with strand here or does project deal with that?


	#Now we have to deal with DBEntries/AssociatedFeatureTypes

	if($conf->param('associations')){
	  $feat->get_all_DBEntries;#this will simply store this in $feat->{'dbentries'}

	  #This is currently only done by default for ExternalFeatures
	  #So would need to explicitly store these for others
	  $feat->associated_feature_types;#This will store in $feat->{'associated_feature_types'}
	}

	$feat->{'dbID'} = undef;
	$feat->{'adaptor'} = undef;
	push @new_feats, $feat;
	$stored_cnt ++;
  }

  #No need to reset dnadb as we have attached the correct one in the feature slice
  #Will this store DBEntries and associated_feature_types
  #This is a balance between calling the store accessor repeatedly here
  #and testing the internal caches for all feature stores
  #Let's do it here for now.

  if (@new_feats){

	foreach my $feat(@new_feats){

	  #Associated ftypes will bestored by default for ExternalFeatures
	  ($feat) = @{$set_feat_adaptor->store($feat)};

	  if($conf->param('associations')){

		#Now restore all DBEntries
		
		foreach my $dbentry(@{$feat->get_all_DBEntries}){
		  
		 
		  #This fails if edb.version is NULL!

		  $dbentry_adaptor->store($dbentry, $feat->dbID, $feat_class);#Do we need to ignore release here?
		  #But this should have the same info as the stored DBEntry
		  #Is this because db_release is not being returned?


		}
	  }
	}
  }

  #test just to be sure
 
  my @tmp = @{$fset->get_Features_by_Slice($new_slice)};
		   

  $logger->info("Projected ${stored_cnt}/${num_old_feats} features from ".$slice->name.' to new assembly '.
				$conf->param('new_assembly')."\n");

  $logger->info("$no_projection_cnt have no mapping\n");
  $logger->info("$multi_segment_cnt mapped to multiple regions\n");
  $logger->info("$wrong_length_cnt $length_txt\n");
  $logger->info("$failed_cnt features failed to project to the new assembly\n");
  $logger->info('Retrieved '.scalar(@tmp).' features from new slice '.$new_slice->name()."\n");
  
  $total_no_proj       += $no_projection_cnt;
  $total_multi_segment += $multi_segment_cnt;
  $total_failed        += $failed_cnt;
  $total_wrong_length  += $wrong_length_cnt;
  $total_stored        += $stored_cnt;

}

$logger->info('Total '.$fset->name.' features projected to new assembly '.$conf->param('new_assembly').":\t${total_stored}/$total_old_feats\n");
$logger->info("Total $total_wrong_length $length_txt\n");
$logger->info("Total multi region mappings:\t $total_multi_segment\n");
$logger->info("Total with no mapping $total_no_proj\n");
$logger->info("Total $total_failed features failed to project to the new assembly\n");
# finish logfile
$logger->finish_log;

1;
