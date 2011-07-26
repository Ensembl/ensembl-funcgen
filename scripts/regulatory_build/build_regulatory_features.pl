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

build_regulatory_features.pl -- builds features for the "Ensembl 
Regulatory Build", the moral equivalent of the gene build

=head1 SYNOPSIS

run_build_regulatory_features.pl -host host -user user -pass password 
    -outdir output_directory -focus feature_setA,feature_setB  
    -target feature_setC,feature_setD -seq_region_name string


Options:

  Mandatory
    -host|h            Host for eFG DB
    -user|u            User for eFG DB
    -pass|p            Password for eFG DB
    -dbname|d          Name of eFG DB
    -data_version|v    Version of data in eFG DB (e.g. 51_36m)
    -outdir|o          Name of outputut directory
    
    -focus|f           Focus features
    -attrib|a          Attribute features
    OR
    -use_tracking_db   Use tracking DB flag, gets input sets from tracking DB

  Optional
    -port              Port for eFG DB, default is 3306
    -species           Latin species name e.g. homo_sapiens

    -slices|s 
    -gene_signature 
    -cell_type_projection      This builds on core regions for all cell types, regardless of
                               whether a focus feature is present for that cell_type.
    -write_features|w
    -rollback
    -dump_annotated_features
    -dump_regulatory_features
    -logfile                  Defines path to log file
    -tee                       Outputs lof to screen as well as log file
    -help              Print this message and exits



=head1 DESCRIPTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! STRATEGY CHANGED; DOCUMENTATION NEEDS TO BE UPDATED !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

This script is the core to compute the regulatory features for the 
"Ensembl Regulatory Build". A regulatory feature consists of

 a) all features that directly overlap with a focus feature, and
 b) features that are contained in focus-feature-overlaping features

The following figure gives examples.

      |------|  |-------- F1 --------|   |----|

      |---------------------------------------|

  |--X--|  |------| |------|   |---------|  |--X---|

      |============= RegFeature ==============|

                                            |- F2 -|

      |=============== RegFeature =================|

[more documentation to be added]

=cut


#To do

# 1 DONE There are warns and STDERR prints, but we don't specify an -e or a -o file when bsubing
#   This has the knock on effect of not knowing what has failed, as we're arbitrarily using the
#   job index to select a chromosome from an array.  As this does not necessarily map to the job name it's 
#   Could also separate the STDERR and STDOUT, not-essential.

# 2 Need to run using non-ref toplevel, currently just using toplevel

# 3 Remove LSBJOBINDEX stuff? Easier to ID output files by slice

# 4 Dump annotated_features using seq_name not seq_region_id

# 5 DONE Implement Helper::rollback_FeatureSet

# 6 Need to check running jobs before resubmitting, otherwise we get duplicates at the end of a chr.
#   As the second job deletes the early written features from the first job

# SG: Need to double-check that feature_set and data_set are stored and have status "displayable" as well as 
#     having supporting sets correctly associated with the data_set. In v52 there was still a problem.

# 7 Add stats on build: mean, median, mode of the core and whole rf regions? Or just do as separate process, so we can rerun at anytime across all chrs?

# 8 Update RegulatoryRegion analysis text, source from analysis.txt?

# 9 DONE Now uses get_seq_region_id_by_Slice rather than unsafe schema_build

# 10 CellType RFs are built on the scaffold even if the only focus present on that cell line exceeds the max length and is demoted to an attr? This would require remembering the cell type of the demoted focus features across all RF scaffolds which is spans.

# 11 DONE Now tests before setting attr end in update_focus

# 12 DONE Fixed AnnotatedFeatures dump bug

# 13 Populate DataSet supporting sets!

# 14 Having absent RFs where there is a core RF may obfuscate the data
#    Is it not better to have a redundant core RF here?

# 15 Fix update_attributes! Is currently not setting attr start/end correctly

# 16 DONE Implement range registry to avoid overlap/encapsulated calc errors

# 17 DONE Now tests bound vs attribute start/end when storing

# 18 ONGOING Improve debug output

# 19 DONE MT is always getting submitted?

# 20 Remove dumps completely and slurp straight into perl, max is about 20 MB of data?

# 21 Check 'removed' focus features are actually being added as attributes?

# 22 DONE Log to file and implement tee

# 23 Look at feature count validation, currently warning instead of warning

# 24 Handle no available features better, die for chromosomes, exit for non chromosomes?

# 25 Fix archiving see ~ line 491. dset can't be true due to previous throw

# 26 meta keys still not being stored? Where should this go? Here, define_and_validate_set or DataSetAdaptor?
#    Now done in DataSetAdaptor, but we currently can't set focus string as this doesn't know about them
#    Set focus/attribute status for feature sets? Overkill? Move back here?

# 27 Build on core regions projected to cell lines which do not have core data.

# 28 Add attribute method

# 29 Don't build on cell type with only spare focus features and no other attrs i.e. only DNase
#    Min of N different focus features if no attrs? Just include in core set.

# 30 Archive can get very messy if it fails half way through as it is non-recoverable i.e.
#    Next run will not try and archive and will treat everything as current.


use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Funcgen::Utils::Helper;
use Bio::EnsEMBL::Utils::Exception qw(verbose warning info);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(get_current_regulatory_input_names open_file strip_param_flags strip_param_args generate_slices_from_names);
use Bio::EnsEMBL::Funcgen::RegulatoryFeature;
use Bio::EnsEMBL::Funcgen::FeatureSet;
local $|=1;

my ($pass,$port,$host,$user,$dbname,
	$dnadb_pass,$dnadb_port,$dnadb_host,$dnadb_user,$dnadb_name,
    $outdir, $write_features, $cdb,
    $no_dump_annotated_features,$dump_regulatory_features, $include_mt,
    $clobber,$focus_max_length, @slices, @skip_slices, $non_ref,
    $focus_extend, $dump,$gene_signature, $ctype_projection, $use_tracking,
    $debug_start, $debug_end, $version, @focus_names, @attr_names, $local);



#Declare to avoid only used once warnings;
$main::_tee      = undef;
$main::_log_file = undef;
my @tmp_args = @ARGV;

#remove these defaults?
$host = $ENV{EFG_HOST};
$port = $ENV{EFG_PORT};
$user = $ENV{EFG_WRITE_USER};
$dbname = $ENV{EFG_DBNAME};
#$species = $ENV{SPECIES};

GetOptions (
            "pass|p=s"       => \$pass,
            "port=s"         => \$port,
            "host|h=s"       => \$host,
            "user|u=s"       => \$user,
            "dbname|d=s"     => \$dbname,
			"dnadb_pass|p=s"       => \$dnadb_pass,
            "dnadb_port=s"   => \$dnadb_port,
            "dnadb_host|h=s" => \$dnadb_host,
            "dnadb_user|u=s" => \$dnadb_user,
            "dnadb_name|d=s" => \$dnadb_name,
		       
			"outdir|o=s"     => \$outdir,
            #"do_intersect|i=s" => \$do_intersect,
            "write_features|w" => \$write_features,
            "no_dump_annotated_features"  => \$no_dump_annotated_features,
            "dump_regulatory_features"  => \$dump_regulatory_features,
			"rollback" => \$clobber,
			"focus_max_length=i" => \$focus_max_length,
            "focus_extend=i" => \$focus_extend,
			"dump" => \$dump,
            "gene_signature" => \$gene_signature,
			"debug_start=i" => \$debug_start,
            "debug_end=i" => \$debug_end, #Defaults to slice end if not set

			'cell_type_projection'=> \$ctype_projection,
			'skip_slices=s{,}'    => \@skip_slices,
			'slices=s{,}'         => \@slices,
			'focus_sets=s{,}'     => \@focus_names, 
			'attribute_sets=s{,}' => \@attr_names,
			'use_tracking_db'     => \$use_tracking,
			'non_ref'             => \$non_ref,
			'include_mt'          => \$include_mt,
			'version=s'           => \$version,
			'local'               => \$local,
			'tee'                 => \$main::_tee,
			'logfile=s'             => \$main::_log_file,
	
			'help|?'         => sub { pos2usage(-exitval => 0, -message => "Params are:\t@tmp_args"); }
		   ) or pod2usage( -exitval => 1);


my $helper = new Bio::EnsEMBL::Funcgen::Utils::Helper();
$helper->log("build_regulatory_features.pl @tmp_args");

die('focus_extend may not safe, check update_attributes') if $focus_extend;


#Allow comma separated quoted names containing spaces
@focus_names = split(/,/o, join(',',@focus_names));
@attr_names  = split(/,/o, join(',',@attr_names));

### defaults ###
$port ||= 3306;

$dnadb_port ||= $port;
$dnadb_host ||= $host;
$dnadb_user ||= $user;
$dnadb_pass ||= $pass;
						
### check options ###

die("Must specify mandatory -host\n") if ! defined $host;
die("Must specify mandatory -user\n") if ! defined $user;
die("Must specify mandatory -pass\n") if ! defined $pass;
die("Must specify mandatory -dbname\n") if ! defined $dbname;
die("Must specify mandatory output directory (-outdir).\n")       if ! $outdir;



$focus_max_length = 2000 if (! defined $focus_max_length);
$focus_extend = 2000 if (! defined $focus_extend);

#die("No output directory specified! Use -o option.") if (!$outdir);
if (defined $outdir && ! -d $outdir) {
  system("mkdir -p $outdir");
}

$outdir =~ s/\/$//o;



#use Bio::EnsEMBL::Funcgen::Utils::RegulatoryBuild qw(is_overlap);
#NJ what was this being used for and can we use the RangeRegistry?


# use ensembldb as we may want to use an old version?
#NJ Default should be staging, but add params for overriding
if ($dnadb_name){
  $cdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
	(
	 -host => $dnadb_host,
	 -port => $dnadb_port,
	 -user => $dnadb_user,
	 -pass => $dnadb_pass,
	 -dbname => $dnadb_name,
	 -group   => 'core',
	);
}

if($gene_signature){
  warn "You need to check that the schema_build of your efg and dnadb match "
}


my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
  (
   -host   => $host,
   -user   => $user,
   -dbname => $dbname,
   -pass   => $pass,
   -port   => $port,
   -dnadb  => $cdb,
   -group  => 'funcgen',#Should be set as default in adaptor new method
  );


#Test connections
$db->dbc->db_handle;
$db->dnadb->dbc->db_handle;

#Test species
my $species = $db->species;

if(! defined $species){
  die("Could not get a valid species from $dbname, please check the meta species.production_name");
}


if(! $use_tracking){
  die("Must specify mandatory focus sets (-focus_sets).\n")         if ! @focus_names;
  die("Must specify mandatory attribute sets (-attribute_sets).\n") if ! @attr_names;
}
elsif(@focus_names || @attr_names){
  die('Please specify either -use_tracking or -focus_sets and -attribute_sets');
}
else{#Use Tracking

  #Change this to Funcgen when we implement full funcgen schema

  my $tdb = Bio::EnsEMBL::DBSQL::DBConnection->new
	(
	 -host => 'ens-genomics1',
	 #-port => $dnadb_port,
	 -user => 'ensro',
	 #-pass => $dnadb_pass,
	 -dbname => 'efg_data_tracking',
	 #-group   => 'core',
	);

  @focus_names = get_current_regulatory_input_names($tdb, $db, 1);#Focus flag
  @attr_names = get_current_regulatory_input_names($tdb, $db);
}


my $fsa = $db->get_FeatureSetAdaptor();
my $dsa = $db->get_DataSetAdaptor();
my $fta = $db->get_FeatureTypeAdaptor();
my $cta = $db->get_CellTypeAdaptor();
my $afa = $db->get_AnnotatedFeatureAdaptor();
my $rfa = $db->get_RegulatoryFeatureAdaptor();
my $efa = $db->get_ExternalFeatureAdaptor();
my $sa  = $db->dnadb->get_SliceAdaptor();
my $aa  = $db->get_AnalysisAdaptor();
my $ga  = $db->dnadb->get_GeneAdaptor();



#Is this true across all supporting and reg sets?
#Do we want to handle the rsets here?
my ($dset_states, $rset_states, $fset_states) = $helper->get_regbuild_set_states($db);




### Validate Focus/Attribute FeatureSets

# parse focus and attribute sets and check that they exist
my (%focus_fsets, %attrib_fsets);

foreach my $fname(@focus_names){
  my $fset = $fsa->fetch_by_name($fname); 
  warn "Fetching focus set:\t$fname\n" if $debug_start;
  die("Focus set $fname does not exist in the DB") if (! defined $fset); 
  $focus_fsets{$fset->dbID} = $fset; 
}
  
foreach my $aname(@attr_names){ 
  my $fset = $fsa->fetch_by_name($aname);
  warn "Fetching attribute set:\t$aname\n" if $debug_start;
  die("Attribute set $aname does not exist in the DB") if (! defined $fset); 
  $attrib_fsets{$fset->dbID()} = $fset; 
}



#This should all be done before the storing(archiving) of the new set
#And before the jobs are submitted to the farm
#Define CellTypes sets

my (%cell_types, %fset_feature_types, %fset_cell_types, %feature_types);

foreach my $fset (values(%focus_fsets)) {
  $cell_types{$fset->cell_type->name} = $fset->cell_type;
  $fset_cell_types{$fset->dbID} = $fset->cell_type->name;
}

#Add special MultiCell ctype
my $mc_ctype = $cta->fetch_by_name('MultiCell');
die ('Your DB does not contain the required MultiCell CellType') if ! $mc_ctype;
$cell_types{MultiCell} = $mc_ctype;

#Validate we don't have any attr set cell types which are 
#not represented in the focus sets

foreach my $fset (values %attrib_fsets) {

  if (! exists $cell_types{$fset->cell_type->name}) {
	
	if(! $ctype_projection){
	  die("You have defined an attribute set whose CellType is not represented in the focus sets:\t".$fset->name);
	}
	$cell_types{$fset->cell_type->name} = $fset->cell_type;
  }

  $fset_cell_types{$fset->dbID} = $fset->cell_type->name;
}



### Archive/Fetch/Store RegulatoryFeature sets.
### CHECK ARCHIVED VERSIONS
#so we don't clobber it before we have a chance to stable_id map
#Or should we remove this facility of host more than one reg build and map between DBs?
#Move this to EFGUtils/Healthchecker for reuse in update_DB_for_release?
#merge with get_regulatory_FeatureSets

if(! $version){
  die('To properly archive a previous regulatory build, you must provide the version number of the new build');
}

my $old_version = $version - 1;
my $mc = $db->get_MetaContainer;
my $some_old_not_archived = 0;
my $some_old_archived = 0;
my @meta_keys = ('regbuild.%s.feature_set_ids', 'regbuild.initial_release_date', 'regbuild.%s.feature_type_ids', 'regbuild.%s.focus_feature_set_ids', 'regbuild.last_annotation_update');

### Validate regbuild version
my ($current_version) = @{$mc->list_value_by_key('regbuild.version')};

if((! defined $current_version) && 
   ($version != 1)){
  die('Could not find regbuild.version meta entry, please correct manually');
}

if(($current_version != $old_version) &&
   ($current_version != $version)){
  die("The regbuild version specified($version) is not the next or current version according to meta regbuild.version($current_version). Please correct manually\n");
}


### ARCHIVE PREVIOUS META/FEATURE/DATASETS
#Move this to EFGUtils::archive_RegulatoryFeatureSet?

if($current_version == $version){#Only archive if version is different to current version
  #Do we need to do some rollback here or maybe in build_reg_features

  if(! $clobber){
	die("A RegulatoryFeatures set already exists with version $version\nPlease update your -version or specify -rollback to overwrite\n");
  }
}
else{


  #This is failing to rollback regbuild.NPC.feature_set_ids?
  #But only happens once as next run will be caught above
  #Is this due to a fail half way through a archive?
  #Need to make this recoverable or stop here is there are any partial archive sets found
  #do we can stop and correct before any more mess is made.
  

  foreach my $ctype(keys %cell_types){

	foreach my $mkey(@meta_keys){
	  #Check both the current and the old versions just in csae we forgot to update the version when running
	  
	  my $meta_key = sprintf($mkey, $ctype);
	  
	  my ($mvalue) = @{$mc->list_value_by_key("${meta_key}_v${version}")};
	  
	  if($mvalue){
		die("Found meta entry for:\t${meta_key}_v${version}\t$mvalue\n".
			"It appears that the regulatory build version you have specified($version) has already been at least partially archived. Maybe you want to set the version param = ".($current_version+1)."?\nPlease correct manually.\n");
	  }
	  
	  if($version != 1){
		($mvalue) = @{$mc->list_value_by_key("${meta_key}_v${old_version}")};
		
		if(! $mvalue){
		  $some_old_not_archived = 1;
		  $helper->log("Found no meta entry for:\t${meta_key}_v${old_version}");
		}
		else{
		  $helper->log("Found meta entry for:\t${meta_key}_v${old_version}\t$mvalue");
		  $some_old_archived = 1;
		}
	  }
	  else{
		$some_old_archived = 1;
	  }
	}

	
	if($some_old_archived && $some_old_not_archived){
	  die("Archiving of the old RegulatoryFeature(v${old_version}) set has not be completed correctly.\n".
		  "Please correct manually.\n");
	}
	else{ #test old and new feature/data set here for sanity
	  #So we aren't depending solely on the meta keys

	  my $set_name = "RegulatoryFeatures:$ctype";
	  my $archd_set_name = "${set_name}_v${version}";
	  my $old_archd_set_name = "${set_name}_v${old_version}";
	  my $sql;
	  my $fset = $fsa->fetch_by_name($archd_set_name);
	  
	  if($fset){
		die("It appears that the new $archd_set_name FeatureSet has already been archived, but no meta entries are present. Please correct manually.\n")
	  }
	  
	  my $dset = $dsa->fetch_by_name($archd_set_name);
	  
	  if($dset){
		die("It appears that the new $archd_set_name DataSet has already been archived, but no meta entries are present. Please correct manually.\n")
	  }
	  
	  
	  if($version != 1){
		$fset = $fsa->fetch_by_name($old_archd_set_name);
		$dset = $dsa->fetch_by_name($old_archd_set_name);
	  }
	  
	  if($some_old_not_archived){#all old meta not archived
		
		#Test old sets
		if(defined $fset){
		  die("It appears that the old $old_archd_set_name FeatureSet has already been archived, but no meta entries are present. Please correct manually.\n")
		}
		
		if(defined $dset){
		  die("It appears that the old $old_archd_set_name DataSet has already been archived, but no meta entries are present. Please correct manually.\n")
		}
		
		if($version != 1){
		  $helper->log("Archiving old $set_name set");
		  
		  #Sanity check that the sets exists
		  $fset = $fsa->fetch_by_name($set_name);
		  $dset = $dsa->fetch_by_name($set_name);
		  

		  
		  
		  if(! defined $fset){
			warn("You are trying to archive $set_name FeatureSet (current version = $old_version) which does not exist and has not been archived");
		  }
		  
		  if(! defined $dset){
			warn("You are trying to archive $set_name DataSet (current version = $version) which does not exist and has not been archived");
		  }
		  
		  if($fset){
			#We shouldn't need to track this process with a status
			#As we have the exhaustive tests above.
			
			#Rollback all other older feature_sets? Using rollback_FeatureSet
			print "Do we need to rollback_FeatureSet for older versions?\n";
			
			#rename data_set
			$sql = "update data_set set name='${old_archd_set_name}' where name='$set_name'";
			$db->dbc->db_handle->do($sql);
			$dset->adaptor->revoke_status('MART_DISPLAYABLE', $dset);
			$dset->adaptor->revoke_status('DISPLAYABLE', $dset);
			#Keep DAS_DISPLAYABLE status?
			
			#rename feature_set
			$sql = "update feature_set set name='${old_archd_set_name}' where name='$set_name'";
			$db->dbc->db_handle->do($sql);
			
			#validate update and revoke states
			$fset = $fsa->fetch_by_name($old_archd_set_name);
			
			if (! $fset) {
			  die("Failed to create archive FeatureSet:\t${old_archd_set_name}\n");
			}
			
			$fset->adaptor->revoke_status('MART_DISPLAYABLE', $fset);
			$fset->adaptor->revoke_status('DISPLAYABLE', $fset);
			#Keep DAS_DISPLAYABLE status?
			
			#rename meta keys
			foreach my $mkey (@meta_keys) {
			  my $meta_key = sprintf($mkey, $ctype);
			  $sql = "update meta set meta_key='${meta_key}_v${old_version}' where meta_key='${meta_key}'";
			  $db->dbc->db_handle->do($sql);
			}
		  }
		}
		
		#Finally update regbuild.version
		$sql = "update meta set meta_value='${version}' where meta_key='regbuild.version'";
		$db->dbc->db_handle->do($sql);
	  } 
	  else {					# old has been archived
		#Will never happen for version == 1
	  
		if (!$fset) {
		  die("It appears that all the old RegulatoryBuild meta entries have been archived, but archived FeatureSet $old_archd_set_name is not present.\nPlease correct manually.\n")
		}
	  
		if (! $dset) {
		  die("It appears that all the old RegulatoryBuild meta entries have been archived, but archived DataSet $old_archd_set_name is not present.\nPlease correct manually.\n")
		}
	  
		$helper->log("Old $set_name set has previously been archived");
	  }
	}
  }
}


# make sure that attribute sets also contain focus sets 
# and ctype_fsets contain core fsets
my %ctype_fsets;

for my $fset_id(keys %focus_fsets){
  $attrib_fsets{$fset_id}        = $focus_fsets{$fset_id};
  push @{$ctype_fsets{MultiCell}}, $focus_fsets{$fset_id};
}


#Set up some more caches for all of the fsets

foreach my $fset(values %attrib_fsets){
  my $ft_name = $fset->feature_type->name;

  #if(exists $ftype_pwm_names{$species}{$ft_name}){
  #	$pwm_ftype_names{$ftype_pwm_names{$species}{$ft_name}} = $ft_name;
  #  }

  $fset_feature_types{$fset->dbID}          = $ft_name;
  $feature_types{$ft_name}                  = $fset->feature_type;
  push @{$ctype_fsets{$fset->cell_type->name}}, $fset;
}


#Set name sorted fset arrays by ctype for bin strings
foreach my $ctype_name(keys %ctype_fsets){
  @{$ctype_fsets{$ctype_name}} = sort {$a->name cmp $b->name} @{$ctype_fsets{$ctype_name}};
};



### Print some details
foreach my $ctype(keys %ctype_fsets){
  next if $ctype eq 'MultiCell';

  my $fset_txt = "\n\n# $ctype Focus Sets:\n";
  my $attr_txt = "\n# $ctype Attribute Sets:\n";

  foreach my $fset(@{$ctype_fsets{$ctype}}){
	
	if($focus_fsets{$fset->dbID}){
	  $fset_txt .= $fset->name.'('.$fset->dbID.')  ';
	}
	else{
	  $attr_txt .= $fset->name.'('.$fset->dbID.')  ';
	}
  }
  $helper->log($fset_txt);
  $helper->log($attr_txt);

}


### Check whether analysis is already stored
#TO DO Update the description text here? Use flat file import?

my $analysis = Bio::EnsEMBL::Analysis->new
  (
   -logic_name      => 'RegulatoryRegion',
   -db              => 'NULL',
   -db_version      => 'NULL',
   -db_file         => 'NULL',
   -program         => 'NULL',
   -program_version => 'NULL',
   -program_file    => 'NULL',
   -gff_source      => 'NULL',
   -gff_feature     => 'NULL',
   -module          => 'NULL',
   -module_version  => 'NULL',
   -parameters      => 'NULL',
   -created         => 'NULL',
   -description     => 'Union of focus features, features overlapping focus features,'.
   ' and features that are contained within those',
   -display_label   => 'RegulatoryRegion',
   -displayable     => 1,
  );

my $logic_name = $analysis->logic_name();
my $ana = $aa->fetch_by_logic_name($logic_name);

if ( ! defined $ana ) {			# NEW
  $aa->store($analysis);
    
} elsif ( $ana->compare($analysis) ) { # EXISTS, but with different options

  ### analysis compare
  # returns  1 if this analysis is special case of given analysis
  # returns  0 if they are equal
  # returns -1 if they are completely different
    
  die('Analysis with logic name \''.$logic_name.'\' already exists, but '.
	  "has different options! Use different logic_name for you analysis '$logic_name'!");

  #$self->efg_analysis->dbID($analysis->dbID);
  #$self->efg_analysis->adaptor($self->efgdb->get_AnalysisAdaptor);
  #$aa->update($self->efg_analysis);
 
}

$analysis = $aa->fetch_by_logic_name($logic_name);

# Do this here before bsubing to avoid race condition
# This will rollback based on slices
# However, we may get partial rollback and redefining of supporting sets!!
# Can only warn about this?


### Define slices and run or bsub

#We need to define rollback slices here
#As rollback will fail for full delete otherwise?

push @skip_slices, 'MT' if ! $include_mt;#Always skip MT by default
#Is this inc dups here 
@slices = @{&generate_slices_from_names($sa, \@slices, \@skip_slices, 'toplevel', $non_ref, 'inc_dups')};


#Hack to update dset stored without reg string
#foreach my $ctype(keys %cell_types){
#  warn $ctype;
#  my $dset_name = ($ctype eq 'core') ? 'RegulatoryFeatures' : "RegulatoryFeatures:{$ctype}";
#  my $dset = $db->get_DataSetAdaptor->fetch_by_name($dset_name);
#  warn "got $dset";
#}
#exit;

#Doing this here for all the slices before we batch will be faster
my $rfsets = &get_regulatory_FeatureSets($analysis, \%cell_types);


#This will get rolled back again unnecessarily for each slice job, 
#but it should still be much quicker
#Well, it will be when rollback_FeatureSet supports multiple slices
my $sr_name;

if (! $local){ #BSUB!
  
  @tmp_args = @{&strip_param_flags(\@tmp_args, ('local', 'non_ref', 'include_mt', 'tee'))};
  @tmp_args = @{strip_param_args(\@tmp_args, ('slices', 'skip_slices'))};
  (my $bsubhost = $host) =~ s/-/_/g;


  foreach my $slice(@slices){
	$sr_name = $slice->seq_region_name;
	my $bsub = "bsub -q normal -o $outdir/RegulatoryBuild_$sr_name.out -e $outdir/RegulatoryBuild_$sr_name.err ".
	  "-J 'RegulatoryBuild_${sr_name}_${dbname}' -R 'select[my".$bsubhost."<80] rusage[my".$bsubhost."=10:duration=10]' -R 'select[mem>4000] rusage[mem=4000]' -M 4000000";
	
	my $cmd = "$ENV{EFG_SRC}/scripts/regulatory_build/build_regulatory_features.pl -local -slices $sr_name @tmp_args";
	  
	$helper->log("Submitting job for slice:\t".$sr_name);
	system("$bsub $cmd") && die("Can't submit job array to farm");#Is this failure getting caught properly?

  }
  
  $helper->log("Output can be found here:\t$outdir");
  warn "We need a to put this in a Runnable so we can easily monitor which jobs have failed";
  exit;

}
elsif(scalar(@slices) > 1){
  die("Cannot run in local mode with more than one slice. Please omit -local or specify a single slice name.\n");
}


### Run slice job locally
my $slice = $slices[0];
$sr_name = $slice->seq_region_name;

# get core and fg seq_region_id for slice
my ($core_sr_id, $fg_sr_id);
$core_sr_id = $slice->get_seq_region_id;

# need to build cache first
$afa->build_seq_region_cache();
$fg_sr_id = $afa->get_seq_region_id_by_Slice($slice);
my $slice_name = $slice->name;
die("No eFG seq_region id defined for ".$slice_name." Almost certain there is no data".
	"available on this slice. You might want to run_update_DB_for_release.pl") 
  if ! defined $fg_sr_id;


### Dump AnnotatedFeatures
my $af_file = $outdir.'/annotated_features.'.$slice->seq_region_name.'.dat';

if ((! -f $af_file) ||
	(! $no_dump_annotated_features)){
  &dump_annotated_features();

  if(! -f $af_file){
	warn "No annotated_features present, aborting build\n";
	exit;
  }
} 
else {
  $helper->log("Using previously dumped annotated features:\t$af_file");
}

my $fh = open_file($af_file);


### Parse and build RegulatoryFeatures
if($debug_start && ! $debug_end){ #Set debug_end if not set
  $debug_end = $slice->end;
}

#Some variables which are used in subs
my (%afs, @rf, $ctype, $af_id, $sr_id, $start, $end, $strand, $score, $fset_id, $mf_ids);
my (%feature_count, %seen_af, %removed_ff);#Stats
my $rf_size = -1;

$helper->log("Processing slice:\t".$slice_name);


while (<$fh>) {
  # Read from file and process sequentially in sorted by start, end order. Each 
  # new feature is checked if it overlaps with the preceeding already seen 
  # features. If yes we just carry on with the next one. Otherwise
  next if (/^\#/o);
  chomp; 

  #This assumes the correct seq_region as we dump out separate file
  #+----------------------+---------------+------+------------------+----------------+-------------------+------------+----------------+----------------+
  #| annotated_feature_id | seq_region_id | name | seq_region_start | seq_region_end | seq_region_strand | score      | feature_set_id | mf_ids         |
  #+----------------------+---------------+------+------------------+----------------+-------------------+------------+----------------+----------------+
  #|              1164424 |           210 | 1    |          4841882 |        4842334 |                 0 |   15.83493 |            112 | 8773,8774      |

  ($af_id, $sr_id, undef, $start, $end, $strand, $score, $fset_id, $mf_ids) = split (/\s+/o, $_);
  
  #Skip non debug regions
  #Warning, this may not build accurately!!!

  #This will not run with a subslice properly
  
  


  if ($debug_start) {
	next if ($start < $debug_start);
	last if ($start > $debug_end);

	print "\nNEW ".$attrib_fsets{$fset_id}->cell_type->name.':'.$attrib_fsets{$fset_id}->feature_type->name."\t${start}\t${end}\n";

  }
    
  #132832790-132833336

  #Peak filtering
  #next if (exists $ChIPseq_cutoff{$attrib_fsets{$fset_id}->name()} 
  #         && $score < $ChIPseq_cutoff{$attrib_fsets{$fset_id}->name()});
  
  $helper->log($_, "\t", $attrib_fsets{$fset_id}->name) if ($debug_start);
  my  $length = $end-$start+1;
  $ctype = $fset_cell_types{$fset_id};
 
  # some stats
  $feature_count{fsets}{$fset_id}++;
  $feature_count{ctypes}{$ctype}++; #Need to report these
	

  if (exists $focus_fsets{$fset_id}) { ### Focus feature
	$feature_count{focus_ctypes}{$ctype}++;

	#Where is the overlap check for this?
	#Let's convert to multi-celltype first then fix
	
	if ($length > $focus_max_length) { # Too long, demote to attribute
	  print "Focus feature ($af_id) longer than $focus_max_length; ",
		"added to attribute feature list\n" if ($debug_start);

	  #Need to build here but not modify core region?
	  #This will become redundant if we all ways do a cell _type projection build

	  &add_feature();

	  #We are not adding this as an attribute feature
	  #Just removing it!!
	  
	  
	  $removed_ff{$fset_id}{$af_id} = 1;
	  $removed_ff{$ctype}++;
	  # need to check overlaps with current focus_features
	  next;
	}
	

		
	# focus feature overlaps w/  rf core (0° of separation)
	if (@rf && 
		($start <= $rf[$rf_size]->{focus_end})){

	  my $update_attrs = 0;

	  if (! exists ${$rf[$rf_size]{fsets}}{$ctype}){
		$update_attrs = 1;
	  }

	  &update_focus();
	  &update_attributes() if $update_attrs;
	} 
	else {					# open new regulatory feature
	  &add_focus();
	  &update_attributes();
	  #How is this capturing further down stream fully contained features???
	  #Are these picked up in the %afs list?
	}
   	&add_feature();
  } 
  else { ### Attribute feature
	
	#Only update if we have seen a focus feature of this ctype
	if( @rf && 
		((exists ${$rf[$rf_size]{fsets}}{$ctype}) || $ctype_projection)){

	  if ($start <= $rf[$rf_size]{focus_end}){
		print "\tAdding overlapping $ctype feature(".$af_id.")\n" if ($debug_start);
	  
		# add annot. feature id to reg. feature
		$rf[$rf_size]{annotated}{$ctype}{$af_id} = $fset_id;
		$rf[$rf_size]{fsets}{$ctype}{$fset_id}++;
		$seen_af{$fset_id}{$af_id} = 1;
		
		# update end of regulatory feature
		#warn "$ctype $end > ".$rf[$rf_size]{attribute_end}{$ctype};

		if ($end > $rf[$rf_size]{attribute_end}{$ctype}){
		  $rf[$rf_size]{attribute_end}{$ctype} = $end ;
		  print "\tResetting $ctype attr end to ".$end."\n" if $debug_start;
		}
		#else this should never happen as the dumps should be sorted by start/end?

	  } 
	  elsif( ($end <= $rf[$rf_size]{attribute_end}{$ctype}) &&
			 ($start <= ($rf[$rf_size]{focus_end} + $focus_extend)) ){	
		
		#Do we not need to account for unset attr ends in here too
		#What about completely abset attrs
		#Need to set bound_end/start in create_regualtory_features?
	
		print "Adding fully contained $ctype feature(".$af_id.")\n" if ($debug_start);
		
		# add annot. feature id to reg. feature
		$rf[$rf_size]{annotated}{$ctype}{$af_id} = $fset_id;
		$rf[$rf_size]{fsets}{$ctype}{$fset_id}++;
		$seen_af{$fset_id}{$af_id} = 1;
	  }
	  else{
		print "Out of range $ctype $start-$end\n" if $debug_start; #level 3?
	  }

	  print "\t$ctype attr start/end now:\t".$rf[$rf_size]{attribute_start}{$ctype}.' - '.$rf[$rf_size]{attribute_end}{$ctype}."\n" if $debug_start;

	  
	}

	&add_feature();
  }
}


$helper->log_header("Skipping storage RegulatoryFeatures") if ! $write_features;

&create_regulatory_features(\@rf);


if ($dump_regulatory_features) {
    
  #Change this to use @reg_feats from create_regulatory_features?

  my $outfile  = $outdir.'/regulatory_features_'.$slice->seq_region_name.'.dat';
  my $out = open_file($outfile, ">");
    
  foreach my $regf(@rf){
	
	printf $out "%d\t%s\t%d\t%d\t%d\t%d\t%s\t%s\n", 
	  $fg_sr_id, $slice->seq_region_name, 
		$regf->{focus_start}, $regf->{focus_end},
		  $regf->{attribute_start}, $regf->{attribute_end}, 
			$regf->{binstring}, join(",", sort {$a<=>$b} keys %{$regf->{annotated}});
  };
  
  $helper->log(sprintf "# Regulatory features written to ".$outfile);
  
}



$helper->log_header("RegulatoryBuild Report");

my $total_feature_count = 0;


#Do this in create_regulatory_features
#my (%rf_count);
#map { 
#  foreach my $k (keys %{$_->{fsets}}) {
#	#print join(" ", $k, $_->{fsets}->{$k}), "\n";
#	$rf_count{$k} += $_->{fsets}->{$k};
#  }
#} @rf;

my $die = 0;

foreach my $ctype(keys %ctype_fsets){
  next if $ctype eq 'MultiCell';

  $helper->log(sprintf "#$ctype FeatureSets %80s\t%8s\t%8s\t%8s",'Total', 'Included', 'Distinct', 'Removed');
  

  foreach my $set(@{$ctype_fsets{$ctype}}){
	my $set_feature_cnt = $feature_count{fsets}{$set->dbID} ||= 0;
	$total_feature_count += $set_feature_cnt;
	
	$helper->log(sprintf "#%80s(%d)\t%8d\t%8d\t%8d\t%8d", $set->name, $set->dbID, 
			   $set_feature_cnt, #Total
				 $feature_count{regfeats}{$set->dbID} || 0, #Included
				 scalar(keys %{$seen_af{$set->dbID}}) || 0, #Distinct
				 scalar(keys %{$removed_ff{$set->dbID}}) || 0); #Removed
  }
  
  my $total = $feature_count{ctypes}->{$ctype} || 0;
  my $focus = $feature_count{focus_ctypes}->{$ctype} || 0;
  $helper->log("");
  $helper->log("# Total $ctype features seen:\t$total");
  $helper->log("# Focus $ctype feature seen:\t$focus");

  if ((($total < $focus) || ($focus == 0)) &&
	 ! $ctype_projection){
	warn("Found inconsistenty between Total($total) and Focus($focus) $ctype features\n");
	$die = 1;
  }

  #This number should be the Focus - those focus which have been removed
  my $f_regfs  = $feature_count{regfeats}{$ctype}           || 0; #Non-projected
  my $t_regfs  = $feature_count{total_regfeats}{$ctype}     || 0; #Inc projected
  my $mfs      = $feature_count{motif_features}{$ctype}     || 0;
  my $p_regfs  = $feature_count{projected_regfeats}{$ctype} || 0;
  my $fp_regfs = $feature_count{failed_projected}{$ctype}   || 0; 

  if($ctype_projection){
	$helper->log("# Projected $ctype RegulatoryFeatures:\t$p_regfs");
	$helper->log("# Failed to project $ctype RegulatoryFeatures:\t$fp_regfs");

	if($t_regfs != ($f_regfs + $p_regfs + $fp_regfs)){
	  warn("WARNING:\tTotal RegFeats is not a sum of Non-project & Projected & Failed");
	  $die = 1;
	}
  }


  $helper->log("# Focus        $ctype RegulatoryFeatures:\t$f_regfs");
  $helper->log("# Total        $ctype RegulatoryFeatures:\t$t_regfs");
  $helper->log("# MotifFeature $ctype inclusions:        \t$mfs");
  $helper->log("\n\n");
  $removed_ff{$ctype} ||= 0;
  my $max_regfs = $focus - $removed_ff{$ctype};

  
  if($t_regfs > $max_regfs){
	#Change this to log if ctype_projection?
	warn("WARNING:\tFound more RegFeats($t_regfs) than included focus features($max_regfs) for $ctype\n");
	$die = 1 if ! $ctype_projection;
  }
  elsif($t_regfs == 0){
	warn("WARNING:\tFound no features for $ctype\n");
	$die = 1;
  }
}


$helper->log(sprintf "# Total number of read features: %10d", $total_feature_count);
$helper->log(sprintf "# Total number of RegulatoryFeatures: %10d", scalar(@rf));

warn("ERROR:\tInconsistencies were found between some of your focus sets and the resultant RegulatoryFeatures, see above.\n") if $die;

###############################################################################
# dump annotated features for given seq_region to file
sub dump_annotated_features{
  #To do, remove this and just slurp directly into perl!


  my @fset_ids = keys %attrib_fsets;
  $helper->log("Dumping ".scalar(@fset_ids)." AnnotatedFeature Sets for slice:\t".$slice_name);

  my $sql = 'SELECT af.annotated_feature_id, af.seq_region_id, sr.name, '.
	                'seq_region_start, seq_region_end, seq_region_strand, score, '.
	                'feature_set_id, mf.mf_ids '.
	          "FROM (select distinct(seq_region_id), name from seq_region where seq_region_id='${fg_sr_id}') sr, ".
				    'annotated_feature af '.
		 'LEFT JOIN (SELECT amf.annotated_feature_id, group_concat(mf1.motif_feature_id) as mf_ids '.
		             'FROM  motif_feature mf1, associated_motif_feature amf '.
			         'WHERE amf.motif_feature_id=mf1.motif_feature_id GROUP by amf.annotated_feature_id) mf '.
				'ON mf.annotated_feature_id=af.annotated_feature_id '.
			 'WHERE sr.seq_region_id=af.seq_region_id '.
			   'AND af.feature_set_id in ('.join(',', @fset_ids).') order by seq_region_start, seq_region_end';
	



  #Do we really need the join to sr, can we not just dump the sr_ids instead of the names?
  #We use the core_sr_id to build the slices anyway
		  
  #Delete previous dumps
  unlink($outdir."/annotated_features.".$slice->seq_region_name.".dat");
  
  #This was appending without removing previous dumps!!!!
  #This will not create a file with no data
  my $command = "echo \"$sql\" ".
	" | mysql -quick -N -h".$host." -P".$port." -u".$user." -p".$pass." ".$dbname.
	  #" | gawk '{if (\$7==".join("||\$7==", @fset_ids).") print }'".
	  #	  " | sort -k3,3 -k4,5n".
        " | gawk '{ print >> \"".$outdir."/annotated_features.".$slice->seq_region_name.".dat\" }'";
    
  $helper->log("# Execute: $command");

  # need to remove existing dump files, since we append to the file
  #system("rm -f $outdir/$dbname.annotated_feature_*.dat") &&
  #        die ("Can't remove files");
    
  #This is not catching mysql error due to the pipe.
  #Could do with counting features and validating vs number of lines

  system($command) && die ("Can't dump data to file in $outdir");
    

  return;
}

###############################################################################

#This should set up a new core RF
#We need to track annotation for all cell_type and then only build
#for those cell type for which we have we have seen core evidence
#We can do this by decoding the fset_ids

sub add_focus{
  print "\nNew $ctype focus:\t$start-$end\n" if ($debug_start);

  push @rf, {
			 focus_start       => $start,
			 focus_end         => $end,
			 attribute_start   => { $ctype => $start },
			 attribute_end     => { $ctype => $end },
			 focus             => { 
								   $ctype    => { $af_id => $fset_id },
								   MultiCell => { $af_id => $fset_id },
								  },
			 fsets             => {
								   MultiCell => {$fset_id => 1},
								   $ctype    => {$fset_id => 1},
								  },
			 motif_feature_ids => {},
			};
  
  #Do this here so we don't have to call $# multiple times.
  $rf_size +=1; 
  

  if(defined $mf_ids){  #Cache MotifFeature ids - don't project
	my @motif_feature_ids = split/,/, $mf_ids;
	
	$rf[$rf_size]->{motif_feature_ids}{$ctype}    ||= [];
	$rf[$rf_size]->{motif_feature_ids}{MultiCell} ||= [];
	#$feature_count{motif_features}{MultiCell}      += scalar(@motif_feature_ids);
	$feature_count{motif_features}{$ctype}         += scalar(@motif_feature_ids);
	
	push @{$rf[$rf_size]->{motif_feature_ids}{$ctype}},    @motif_feature_ids;
	push @{$rf[$rf_size]->{motif_feature_ids}{MultiCell}}, @motif_feature_ids;
  }


  if($ctype_projection){   #Set ctype projected attr start/ends
	foreach my $ct(keys %cell_types){
	  next if $ct eq 'MultiCell';
	  next if $ct eq $ctype;
	  print "Initialising $ct projection attr start/end:\t$start - $end\n" if $debug_start;
	  $rf[$rf_size]{attribute_start}{$ct} = $start;
	  $rf[$rf_size]{attribute_end}{$ct}   = $end;

	  #Set this to 0 such that update_focus will also update the projected regfeats
	  $rf[$rf_size]{fsets}{$ct}{$fset_id} ||= 0;

	  #How are we storing only those with attrs?
	}
  }

  $seen_af{$fset_id}{$af_id} = 1;

  return;
}

sub update_focus{
  print "\tUpdating focus(".join(', ', keys(%{$rf[$rf_size]{fsets}})).") with $ctype\n" if ($debug_start);

  #Set start values if we haven't seen this ctype before
  #update_attributes will be called immediately after this
  if(! exists  ${$rf[$rf_size]{fsets}}{$ctype}){
	#Account for previous focus start
	$rf[$rf_size]{attribute_start}{$ctype} = $rf[$rf_size]{focus_start};
	$rf[$rf_size]{attribute_end}{$ctype}   = ($end > $rf[$rf_size]{focus_end}) ? $end : $rf[$rf_size]{focus_end};
	print "Setting $ctype attr start/end to:\t".$rf[$rf_size]{attribute_start}{$ctype}.' - '.$rf[$rf_size]{attribute_end}{$ctype}."\n" if $debug_start;
  }


  $rf[$rf_size]{focus_end} = $end if ($end > $rf[$rf_size]{focus_end});

 
  #Reset attr end for all ctypes if required
  #This prevents bound_seq_region_end < seq_region_end bug
  #which is not visible in the browser
   
  foreach my $rf_ctype(keys %{$rf[$rf_size]{fsets}}){

	next if $rf_ctype eq 'MultiCell';    #no attrs for multi cell

	if(exists $rf[$rf_size]{attribute_end}{$rf_ctype}){

	  if($end   > $rf[$rf_size]{attribute_end}{$rf_ctype}){
		$rf[$rf_size]{attribute_end}{$rf_ctype}   = $end;
		print "Resetting $rf_ctype attr end to new $ctype focus end:\t$end\n" if $debug_start;
	  }
	}
  }


  if(defined $mf_ids){  #Cache MotifFeature ids - don't project
	my @motif_feature_ids = split/,/, $mf_ids;
	
	$rf[$rf_size]->{motif_feature_ids}{$ctype}    ||= [];
	$rf[$rf_size]->{motif_feature_ids}{MultiCell} ||= [];
	#	$feature_count{motif_features}{MultiCell}      += scalar(@motif_feature_ids);
	$feature_count{motif_features}{$ctype}         += scalar(@motif_feature_ids);
	
	push @{$rf[$rf_size]->{motif_feature_ids}{$ctype}},    @motif_feature_ids;
	push @{$rf[$rf_size]->{motif_feature_ids}{MultiCell}}, @motif_feature_ids;
  }

  $rf[$rf_size]{fsets}{MultiCell}{$fset_id}++;
  $rf[$rf_size]{fsets}{$ctype}{$fset_id}++;
  $rf[$rf_size]{focus}{$ctype}{$af_id} = $fset_id;
  $rf[$rf_size]{focus}{MultiCell}{$af_id} = $fset_id;
  


  
  #Do we need to alter this to be celltype aware?
  $seen_af{$fset_id}{$af_id} = 1;

  #print Dumper @rf if ($debug_start);
  return;
}

sub update_attributes{
  my @cts = ($ctype);
  @cts = keys(%cell_types) if $ctype_projection;
  

  #This needs to be called only when we see a new ctype focus
  #Do we build on a ctype is there are encapsulated focus features
  #but non for a given ctype at this locus?
  
  # look upstream for all features overlapping focus feature


  foreach my $ct(@cts){
	next if $ct eq 'MultiCell';
	print "\tUpdating $ct attributes\n" if ($debug_start);
	my $updated = 0;


	for (my $i=0; $i<=$#{$afs{$ct}}; $i++) {
	  
	  if ($afs{$ct}->[$i]{end} >= $rf[$rf_size]{focus_start}) {
		$updated = 1;

		if(scalar(keys %{$rf[$rf_size]{focus}{$ctype}}) == 0){
		  print "Updating projected feature\n" if $debug_start;
		  exit;
		}


		# FIRST update attribute start of the regulatory feature...
		
		if ($afs{$ct}->[$i]{start} < $rf[$rf_size]{attribute_start}{$ct}) {
		  #This is always the case the first feature
		  #If there is one, it will always be upstream or start at the same location as the new focus
		  #And we process the rest in this block so no real need for this test???

		  print "\tUpdating attribute start to:\t\t(".$attrib_fsets{$afs{$ct}->[$i]{fset_id}}->feature_type->name.":".$afs{$ct}->[$i]{af_id}.")\t".$afs{$ct}->[$i]{start}.' (<'.$rf[$rf_size]{attribute_start}{$ct}.")\n" if $debug_start;

		  $rf[$rf_size]{attribute_start}{$ct} = $afs{$ct}->[$i]{start};
		  $rf[$rf_size]{annotated}{$ct}{$afs{$ct}->[$i]{af_id}} = $afs{$ct}->[$i]{fset_id};
		  $rf[$rf_size]{fsets}{$ct}{$afs{$ct}->[$i]{fset_id}}++;


		}
		
		# THEN remove features that are out of scope
		splice(@{$afs{$ct}}, 0, ($i));
			 
		# NOW add the remaining overlapping attribute features and set att end
		foreach my $af(@{$afs{$ct}}) {
		  
		  if ($af->{end} >= $start ||
			  ( $af->{start} >= ($rf[$rf_size]{focus_start} - $focus_extend) )) {
			
			if ($debug_start){
			  print "\tUpdating overlapping attribute feature:\t(".
				$attrib_fsets{$af->{fset_id}}->feature_type->name.
				  ":".$af->{af_id}.")\t".$af->{start}.' - '.$af->{end}."\n";
			}
			
			$rf[$rf_size]{annotated}{$ct}{$af->{af_id}} = $af->{fset_id};
			$rf[$rf_size]{fsets}{$ct}{$af->{fset_id}}++;
			
			#Reset the attr end
			if ($af->{end}   > $rf[$rf_size]{attribute_end}{$ct}) {
			  print "\tUpdating attribute end to:\t\t".$af->{end}.' (> '.$af->{end}.")\n" if $debug_start;
			  $rf[$rf_size]{attribute_end}{$ct} = $af->{end};
			}
			# could also test for demoted focus features here!
			
			$seen_af{$fset_id}{$af_id} = 1;
		  }
		}
		
		last;
	  }
	}

	print "$ct updated coords:\t\t\t\t\t\t".$rf[$rf_size]{attribute_start}{$ct}.' - '.$rf[$rf_size]{attribute_end}{$ct}."\n" if $debug_start;
	print "No upstream $ct features found\n" if $debug_start && (! $updated);
  }

  return;
}


#Need add attribute method

sub add_feature{
  
  print "Adding $ctype feature($af_id) to attribute feature_list\n" if $debug_start;


  push(@{$afs{$ctype}}, 
	   {
		af_id => $af_id,
		start => $start,
		end => $end,
		strand => $strand,
		score => $score,
		fset_id => $fset_id
	   });

  return;
}

sub build_binstring{
  my ($rf, $ct) = @_;
  
  my $binstring = '';
  
  foreach my $fs (@{$ctype_fsets{$ct}}) {
	$binstring .= 
	  (exists $rf->{fsets}->{$ct}->{$fs->dbID})? 1 : 0;	
  }
  
  
  if ($gene_signature) {
	my $tss = 0;
	my $tts = 0;
	
	# get feature slice +/- 2.5kB
	my $expand = 2500;
  
	my ($feature_slice, $expanded_slice);
	$feature_slice = $sa->fetch_by_seq_region_id($core_sr_id,
												 $rf->{focus_start},
												 $rf->{focus_end});
	$expanded_slice = $feature_slice->expand($expand, $expand);
	
	foreach my $g (@{$ga->fetch_all_by_Slice($expanded_slice)}) {
		
	  $tss = 1 if (($g->strand == 1  && $g->start > 0) ||
				   ($g->strand == -1 && $g->end <= $expanded_slice->length));
	  $tts = 1 if (($g->strand == 1  && $g->end <= $expanded_slice->length) ||
				   ($g->strand == -1 && $g->start > 0));
	  
	  last if ($tss && $tts);
	}
	
	$binstring .=  $tss.$tts;
  }
    
  return $binstring;
}

#This is only ever used once!!
#Remove sub or move to module

sub create_regulatory_features{
  my $rfs = shift;

  my @reg_feats;
  print "\n\nCreating ".scalar(@$rfs)." RegulatoryFeatures\n" if $debug_start;


  foreach my $rf (@{$rfs}) {
	print "\nCreating RegulatoryFeature:\t\t\t ".$rf->{focus_start}.' - '.$rf->{focus_end}."\n" if $debug_start;
	my $attr_cache;


	if ($slice->start != 1 || $slice->strand != 1) {
	  warn("**** SLICE doesn't start at position 1; resetting start of slice to 1 ****");
	  $slice->{'start'} = 1;
	  #This only works as we are working with seq_region_start/end values in the dumps
	  #We should really fetch the slices with the include_duplcates flag i.e. 1 as start for PAR regions
	}

	my $focus_start = $rf->{focus_start};
	my $focus_end   = $rf->{focus_end};

	#Now build each ctype RegFeat

	#order cell type so we deal with MultCell first
	#This way we can be sure to catch & cache all the tfbs
	#for inclusion in the other cell types
	#Remove this now we handle MFs differently.
	#my @ctypes = grep { /[^(MultiCell)]/ } keys %{$rf->{fsets}};
	#unshift @ctypes, 'MultiCell';										
	
	foreach my $ct (keys %{$rf->{fsets}}){#(@ctypes) {
	  my $projected = 0;

	  if($debug_start){

		if($ct eq 'MultiCell'){
		  $rf->{attribute_start}{$ct} = $focus_start;
		  $rf->{attribute_end}{$ct}   = $focus_end;
		}
		print "\nBuilding $ct RegFeat with attribute start/ends:\t".
		  join(' - ', ($rf->{attribute_start}{$ct}, "\t\t", "\t\t", $rf->{attribute_end}{$ct}))."\n";
	  }

	  #Only count this is it is a true regf and not a scaffold projection
	  $feature_count{total_regfeats}{$ct}++;
	  
	  #Recording all fset features here, no split on projected
	  foreach my $fset_id(keys %{$rf->{fsets}{$ct}}){
		$feature_count{regfeats}{$fset_id} += $rf->{fsets}{$ct}{$fset_id};
	  }

	  #Build attr cache dependant on presence of attrs
	  #We may not have focus data here if we are doing a ctype_projection build
	  #attr cache is af_id fset_id pairs

	  if(exists $rf->{annotated}->{$ct} && exists $rf->{focus}->{$ct}){
		print "Found annotated and focus attrs\n" if $debug_start;

		$feature_count{regfeats}{$ct}++;
		$attr_cache = {
					   annotated => [keys %{$rf->{focus}->{$ct}},
									 keys %{$rf->{annotated}->{$ct}}],
					  };
	  }
	  elsif( exists $rf->{focus}->{$ct}){
		print "Found focus attrs only\n" if $debug_start;
		#MultiCell and focus only
		$feature_count{regfeats}{$ct}++;
		$attr_cache = { annotated => [keys %{$rf->{focus}->{$ct}}] };
	  }
	  elsif(exists $rf->{annotated}->{$ct}){
		print "Found projected annotated attrs only\n" if $debug_start;
		$projected = 1;
		$feature_count{projected_regfeats}{$ct}++;
		$attr_cache =  { annotated => [keys %{$rf->{annotated}->{$ct}}] };
	  }
	  else{
		#No attrs - This is a failed projection!
		print "Skipping Failed Projection - no attrs found\n" if $debug_start;
		$feature_count{failed_projected}{$ct}++;
		next; #$ct
	  }


	  #Add MotifFeatures as attrs
	  $attr_cache->{motif} = $rf->{motif_feature_ids}{$ctype};
		  
	  #We now need to clean the fset_id values from the attrs cache
	  #map { $attr_cache->{$_} = undef } keys %$attr_cache; 

	  print "Creating $ct RegulatoryFeature with attribute cache:\n".Data::Dumper::Dumper($attr_cache)."\n" if $debug_start;#level 3?

	  my $reg_feat = Bio::EnsEMBL::Funcgen::RegulatoryFeature->new
		(
		 -slice            => $slice,
		 -start            => $focus_start,
		 -end              => $focus_end,
		 -strand           => 0,
		 -binary_string    => &build_binstring($rf, $ct),
		 -feature_set      => $rfsets->{$ct},
		 -feature_type     => $rfsets->{$ct}->feature_type,
		 -attribute_cache  => $attr_cache,
		 -projected        => $projected,
		);
	  
	  if ($write_features) {	#Store
		  
		($reg_feat) = @{$rfa->store($reg_feat)};
		  
		#Sanity check bound here to make sure we have built correctly
		#Can only test bounds if adaptors are set

		if($ct ne 'MultiCell'){#Only have bound for non-core features

		  #Account for bounds being set by focus features from another cell type
		  #How do we account for long 'removed' features, which are added as attrs but do not contribute to bounds/core


		  #This is fetching all attribute features from DB!
		  

		  if ((($reg_feat->bound_start != $rf->{attribute_start}->{$ct}) && ($reg_feat->bound_start != $reg_feat->start)) ||
			  ($reg_feat->bound_end != $rf->{attribute_end}->{$ct}) && ($reg_feat->bound_end != $reg_feat->end)){
			
			warn "$ct attributes:\n\t".
			  join("\t\n",  (
							 map { $_->cell_type->name.':'.$_->feature_type->name.':'.$_->start.':'.$_->end }
							 @{$reg_feat->regulatory_attributes('annotated')} 
							)
				  )."\t\n".
					join("\t\n",  (
								   map { $_->display_label.':'.$_->start.':'.$_->end }
								   @{$reg_feat->regulatory_attributes('motif')}
								  ))
					  ."\n";

				   

			warn "$ct bound_start ".$reg_feat->bound_start.' != seen attr '.$rf->{attribute_start}->{$ct}."\n";
			warn "$ct bound_end   ".$reg_feat->bound_end .' != seen attr '.$rf->{attribute_end}->{$ct}."\n";
			warn "$ct start end   ".$reg_feat->start.' - '.$reg_feat->end."\n";
			die('Calculated attribute_start/end values do not match bound_start/end/focus values');
		  }
		}
	  }
		
	  push @reg_feats, $reg_feat;
	}
  }

  return \@reg_feats;
}
  


sub get_regulatory_FeatureSets{
  my ($analysis, $ctypes) = @_;
  my %rf_sets;
  my $ftype = $fta->fetch_by_name('RegulatoryFeature');
	
  if (! $ftype) {
	  
	$ftype = Bio::EnsEMBL::Funcgen::FeatureType->new
	  (
	   -name        => 'RegulatoryFeature',
	   -description => 'Ensembl Regulatory Feature',
	   -class       => 'Regulatory Feature',
	  );
	
	($ftype) = @{$fta->store($ftype)} if ($write_features);
  }
  
  my (@dsets, @fsets);
	

  foreach my $ctype (keys %{$ctypes}) {	
	my ($desc, $dlabel, $fset_name);
	
	$fset_name = "RegulatoryFeatures:$ctype";
	$dlabel    = "Reg.Feats $ctype";

	if($ctype eq 'MultiCell'){
	  $desc = 'Generic RegulatoryFeature focus regions';
	}
	else{
	  $desc = "$ctype specific RegulatoryFeatures";
	  
	}

	my $rollback = ($clobber) ? 'supporting_sets' : undef;
	
	$helper->log("Defining FeatureSet:\t$fset_name");

	my $dset = $helper->define_and_validate_sets($fset_name, 
												 $analysis, 
												 $ftype, 
												 $ctypes->{$ctype}, 
												 'regulatory',
												 undef,#append flag
												 $db,
												 $ctype_fsets{$ctype},#supporting sets
												 $desc,#desc
												 'sets',#rollback mode
												 1,#recovery?
												 \@slices, $dlabel);

	#Always overwrite in case we have redefined the sets
	&store_regbuild_meta_strings($dset, 1);
	$rf_sets{$ctype} = $dset->product_FeatureSet;

	#Set states
	#Move to Utils/RegBuilder.pm?
   
	push @dsets, $dset;

	foreach my $fset(@{$ctype_fsets{$ctype}}){
	  my $ss_dset = $dsa->fetch_by_product_FeatureSet($fset);

	  if(! $ss_dset){
		die("Could not find DataSet for FeatureSet:\t".$fset->name);
	  }
	  
	  push @dsets, $ss_dset;
	}

	push @fsets, ($dset->product_FeatureSet, @{$ctype_fsets{$ctype}});	
  }

  #Set states
  #Move to Utils/RegBuilder.pm?
  foreach my $dset(@dsets){
	foreach my $ds_state(@{$dset_states}){
	  $dsa->store_status($ds_state, $dset);
	}
  }
    
  foreach my $fset(@fsets){
	foreach my $fs_state(@{$fset_states}){
	  $fsa->store_status($fs_state, $fset);
	}
  }


  $helper->log("Got RegulatoryFeature sets for CellTypes:\t".join(', ', keys %rf_sets));
  return \%rf_sets;

}


#This could move to DataSetAdaptor and be called from define_and_validate_sets
#If focus set info was stored in data set

sub store_regbuild_meta_strings{
  my ($dset, $overwrite) = @_;

  my $ds_adaptor = $dset->adaptor;
  $ds_adaptor->db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::DataSet', $dset);
  my ($sql, $meta_value, $reg_string, $cmd, $msg);
  my $fset = $dset->product_FeatureSet;

  if(! defined $fset ||
	 $fset->feature_class ne 'regulatory'){
	die('You must provide a DataSet with an associated \'regulatory\' product FeatureSet');
  }

  my @ssets = @{$dset->get_supporting_sets};

  if(! @ssets){
	throw('You must provide a DataSet with associated supporting sets');
  }


  my $ctype = (defined $fset->cell_type) ?  $fset->cell_type->name : 'core';

  ### build and import regbuild strings by feature_set_id and feature_type_id

  my $sth = $ds_adaptor->db->dbc->prepare("select meta_value from meta where meta_key='regbuild.${ctype}.feature_set_ids'");
  $sth->execute();
  ($meta_value) = map { "@$_" } @{$sth->fetchall_arrayref};
  $reg_string = join(',', map {$_->dbID} sort {$a->name cmp $b->name} @ssets);

  my %reg_strings = 
	(
	 "regbuild.${ctype}.feature_set_ids" => join(',', map {
	   $_->dbID} sort {$a->name cmp $b->name
					 } @ssets),
	 "regbuild.${ctype}.feature_type_ids" => join(',', map {
	   $_->feature_type->dbID} sort {$a->name cmp $b->name
								   } @ssets),
	);

  
  my @ffset_ids;

  foreach my $ffset_id(keys %focus_fsets){

	if ($focus_fsets{$ffset_id}->cell_type->name eq $dset->cell_type->name){
	  push @ffset_ids, $ffset_id;
	}
  }

  $reg_strings{"regbuild.${ctype}.focus_feature_set_ids"} = join(', ', @ffset_ids);
  


  foreach my $meta_key(keys %reg_strings){
	my $sth = $ds_adaptor->db->dbc->prepare("select meta_value from meta where meta_key='${meta_key}'");
	$sth->execute();
	($meta_value) = map { "@$_" } @{$sth->fetchall_arrayref};


	#warn "meta value for $meta_key is $meta_value";


	if (! defined $meta_value) {
	  $sql = "insert into meta (meta_key, meta_value) values ('${meta_key}', '$reg_strings{${meta_key}}')";
	
	  eval { $ds_adaptor->db->dbc->do($sql) };
	  die("Couldn't store $meta_key in meta table.\n$@") if $@;
	} 
	else {

	  if($meta_value ne $reg_strings{$meta_key}){
		$msg = "$meta_key already exists and does not match\nOld\t$meta_value\nNew\t$reg_strings{${meta_key}}\nPlease archive previous RegulatoryBuild.\n";
		die $msg if(! $overwrite);
		warn $msg;
	  }
	  else{  
		warn "Found matching $meta_key meta entry, do you need to archive a previous Regulatorybuild?\n" if ! $overwrite;
	  }
	}
  }
  
  return \%reg_strings;
}



1;
