#!/software/bin/perl


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
    

  Optional
    -port              Port for eFG DB, default is 3306
    -species           Latin species name e.g. homo_sapiens

    -slices|s 
    -gene_signature 

    -write_features|w
    -clobber
    -dump_annotated_features
    -dump_regulatory_features
    -tee We are already loggin
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

=head1 LICENSE

  Copyright (c) 1999-2009 The European Bioinformatics Institute and
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

# 16 Implement range registry to avoid overlap/encapsulated calc errors

# 17 DONE Now tests bound vs attribute start/end when storing

# 18 ONGOING Improve debug output

# 19 MT is always getting submitted?

# 20 Remove dumps completely and slurp straight into perl, max is about 20 MB of data?

# 21 Check 'removed' focus features are actually being added as attributes?

# 22 Log to file and implement tee

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Funcgen::Utils::Helper;
use Bio::EnsEMBL::Utils::Exception qw(verbose warning info);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(open_file strip_param_flags strip_param_args generate_slices_from_names);
use Bio::EnsEMBL::Funcgen::RegulatoryFeature;
$|=1;

my ($pass,$port,$host,$user,$dbname, $help,
	$dnadb_pass,$dnadb_port,$dnadb_host,$dnadb_user,$dnadb_name,
    $outdir,$do_intersect,$write_features, $cdb,
    $no_dump_annotated_features,$dump_regulatory_features, $include_mt,
    $clobber,$focus_max_length, @slices, @skip_slices, $non_ref,
    $focus_extend, $dump,$gene_signature,
    $debug_start, $debug_end, $version, @focus_names, @attr_names, $local);

my @tmp_args = @ARGV;

#Need to change all print to use Helper log?
my $helper = new Bio::EnsEMBL::Funcgen::Utils::Helper(no_log => 1); ###Turn log file output off for now

print "build_regulatory_features.pl @tmp_args\n";

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
			#"species=s"      => \$species,
       
			"outdir|o=s"     => \$outdir,
            "do_intersect|i=s" => \$do_intersect,
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

			
			'skip_slices=s{,}'    => \@skip_slices,
			'slices=s{,}'         => \@slices,
			'focus_sets=s{,}'     => \@focus_names, 
			'attribute_sets=s{,}' => \@attr_names,
			'non_ref'             => \$non_ref,
			'include_mt'          => \$include_mt,
			'version=s'           => \$version,
			'local'               => \$local,
	
			'help|?'         => sub { pos2usage(-exitval => 0, -message => "Params are:\t@tmp_args"); }
		   ) or pod2usage( -exitval => 1);

die('focus_extend may not safe, check update_attributes') if $focus_extend;


#Allow comma separated quoted names containing spaces
@focus_names = split(/,/,join(',',@focus_names));
@attr_names = split(/,/,join(',',@attr_names));

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
die("Must specify mandatory focus sets (-focus_sets).\n")         if ! @focus_names;
die("Must specify mandatory attribute sets (-attribute_sets).\n") if ! @attr_names;

$focus_max_length = 2000 if (! defined $focus_max_length);
$focus_extend = 2000 if (! defined $focus_extend);

#die("No output directory specified! Use -o option.") if (!$outdir);
if (defined $outdir && ! -d $outdir) {
  system("mkdir -p $outdir");
}
$outdir =~ s/\/$//;




#use Bio::EnsEMBL::Funcgen::Utils::RegulatoryBuild qw(is_overlap);
#NJ what was this being used for and can we use the RangeRegistry?


# use ensembldb as we may want to use an old version
#NJ Default should be staging, but add params for overriding
if ($dnadb_name){
  $cdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
	(
	 -host => $dnadb_host,
	 -port => $dnadb_port,
	 -user => $dnadb_user,
	 -pass => $dnadb_pass,
	 -dbname => $dnadb_name,
	 #-species => $species, 
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
   #-species => $species, #should be set by species.procution_name
   -pass   => $pass,
   -port   => $port,
   -dnadb  => $cdb,
   -group  => 'funcgen',		#Should be set as default in adaptor new method
  );


my $fsa = $db->get_FeatureSetAdaptor();
my $dsa = $db->get_DataSetAdaptor();
my $fta = $db->get_FeatureTypeAdaptor();
my $afa = $db->get_AnnotatedFeatureAdaptor();
my $rfa = $db->get_RegulatoryFeatureAdaptor();
my $sa = $db->dnadb->get_SliceAdaptor();
my $aa = $db->get_AnalysisAdaptor();
my $ga = $db->dnadb->get_GeneAdaptor();



### Validate Focus/Attribute FeatureSets

# parse focus and attribute sets and check that they exist
my (%focus_fsets, %attrib_fsets);
map { my $fset = $fsa->fetch_by_name($_); 
	  warn "Fetching focus set:\t$_\n" if $debug_start;
      die("Focus set $_ does not exist in the DB") if (! defined $fset); 
      $focus_fsets{$fset->dbID} = $fset; 
	} @focus_names;
#print Dumper %focus_fsets;

map { 
  my $fset = $fsa->fetch_by_name($_);
   warn "Fetching attribute set:\t$_\n" if $debug_start;
  die("Attribute set $_ does not exist in the DB") if (! defined $fset); 
  $attrib_fsets{$fset->dbID()} = $fset; 
} @attr_names;



#This should all be done before the storing(archiving) of the new set
#And before the jobs are submitted to the farm
#Define CellTypes sets

my %cell_types;
my %fset_cell_types;#Do we need thsi, can we not jsut dump feature_type_id?

foreach my $fset (values(%focus_fsets)) {
  $cell_types{$fset->cell_type->name} = $fset->cell_type;
  $fset_cell_types{$fset->dbID} = $fset->cell_type->name;
}

$cell_types{core} = undef;

#Validate we don't have any attr set cell types which are 
#not represented in the focus sets

foreach my $fset (values %attrib_fsets) {

  if (! exists $cell_types{$fset->cell_type->name}) {
	die("You have defined an attribute set whose CellType is not represented in the focus sets:\t".$fset->name);
  }

  $fset_cell_types{$fset->dbID} = $fset->cell_type->name;
}




### Archive/Fetch/Store RegulatoryFeature sets.
### CHECK ARCHIVED VERSIONS
#so we don't clobber it before we have a chance to stable_id map
#Or should we remove this facility of host more than one reg build and map between DBs?
#Move this to EFGUtils/Healthchecker for reuse in update_DB_for_release?

if(! $version){
  die('To properly archive a previous regulatory build, you must provide the version number of the new build');
}

my $old_version = $version - 1;
my $mc = $db->get_MetaContainer;
my $some_old_not_archived = 0;
my $some_old_archived = 0;
my @meta_keys = ('regbuild.%s.feature_set_ids', 'regbuild.initial_release_date', 'regbuild.%s.feature_type_ids', 'regbuild.last_annotation_update');

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

  foreach my $ctype(keys %cell_types){


	for my $mkey(@meta_keys){
	  #Check both the current and the old versions just in csae we forgot to update the version when running
	
	  $mkey = sprintf( $mkey, $ctype);


	  my ($mvalue) = @{$mc->list_value_by_key("${mkey}_v${version}")};
	
	if($mvalue){
	  die("Found meta entry for:\t${mkey}_v${version}\t$mvalue\n".
		  "It appears that the regulatory build version you have specified($version) has already been at least partially archived. Maybe you want to set the version param = ".($current_version+1)."?\nPlease correct manually.\n");
	}
	
	if($version != 1){
	  ($mvalue) = @{$mc->list_value_by_key("${mkey}_v${old_version}")};
	
	  if(! $mvalue){
		$some_old_not_archived = 1;
		#print here as this
		print "Found no meta entry for:\t${mkey}_v${old_version}\n";
	  }
	  else{
		print "Found meta entry for:\t${mkey}_v${old_version}\t$mvalue\n";
		$some_old_archived = 1;
	  }
	}
	else{
	  $some_old_archived = 1;
	}
  }

  
  if($some_old_archived && $some_old_not_archived){
	die("Arching of the old RegulatoryFeature(v${old_version}) set has not be completed correctly.\n".
		"Please correct manually.\n");
  }
  else{ #test old and new feature/data set here for sanity
	#So we aren't depending solely on the meta keys

	my $set_name = ($ctype eq 'core') ? 'RegulatoryFeatures' : "RegulatoryFeatures:$ctype";
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
	  if($fset){
		die("It appears that the old $old_archd_set_name FeatureSet has already been archived, but no meta entries are present. Please correct manually.\n")
	  }
	  
	  if($dset){
		die("It appears that the old $old_archd_set_name DataSet has already been archived, but no meta entries are present. Please correct manually.\n")
	  }
	  
	  
	  if($version != 1){
		print "Achiving old RegulatoryFeatures set\n";
		#We shouldn't need to track this process with a status
		#As we have the exhaustive tests above.
	  
		#Rollback all other older feature_sets? Using rollback_FeatureSet
		print "Do we need to rollback_FeatureSet for older versions?\n";
		
		#rename data_set
		$sql = "update data_set set name='${old_archd_set_name}' where name='RegulatoryFeatures'";
		$db->dbc->db_handle->do($sql);
		$dset->adaptor->revoke_status('MART_DISPLAYABLE', $dset);
		$dset->adaptor->revoke_status('DISPLAYABLE', $dset);
		#Keep DAS_DISPLAYABLE status?

		#rename feature_set
		$sql = "update feature_set set name='${old_archd_set_name}' where name='RegulatoryFeatures'";
		$db->dbc->db_handle->do($sql);
		
		#validate update and revoke states
		my $fset = $fsa->fetch_by_name($old_archd_set_name);
		
		if(! $fset){
		  die("Failed to create archive FeatureSet:\t${old_archd_set_name}\n");
		}
	  
		$fset->adaptor->revoke_status('MART_DISPLAYABLE', $fset);
		$fset->adaptor->revoke_status('DISPLAYABLE', $fset);
		#Keep DAS_DISPLAYABLE status?

		#rename meta keys
		foreach my $mkey(@meta_keys){
		  $mkey = sprintf($mkey, $ctype);
		  $sql = "update meta set meta_key='${mkey}_v${old_version}' where meta_key='${mkey}'";
		  $db->dbc->db_handle->do($sql);
		}
	  
		#Finally update regbuild.version
		$sql = "update meta set meta_value='${version}' where meta_key='regbuild.version'";
		$db->dbc->db_handle->do($sql);
	  }
	}
	else{# old has been archived
	  #Will never happen for version == 1
	  
	  if(!$fset){
		die("It appears that all the old RegulatoryBuild meta entries have been archived, but archived FeatureSet $old_archd_set_name is not present.\nPlease correct manually.\n")
	  }
	  
	  if(! $dset){
		die("It appears that all the old RegulatoryBuild meta entries have been archived, but archived DataSet $old_archd_set_name is not present.\nPlease correct manually.\n")
	  }
	  
	  print "Old RegulatoryFeature set has previously been archived\n";
	}
  }
  }
}


# make sure that attribute sets also contain focus sets 
# and ctype_fsets contain core fsets
my %ctype_fsets;
map { 
  $attrib_fsets{$_} = $focus_fsets{$_};
  push @{$ctype_fsets{core}}, $focus_fsets{$_};
} keys %focus_fsets;



#Set name sorted fset arrays by ctype for bin strings
map { push @{$ctype_fsets{$_->cell_type->name}}, $_ } values %attrib_fsets;
map { @{$ctype_fsets{$_}} = sort {$a->name cmp $b->name} @{$ctype_fsets{$_}} } keys %ctype_fsets;

### Print some details
foreach my $ctype(keys %ctype_fsets){
  next if $ctype eq 'core';

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

  print $fset_txt.$attr_txt;

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
@slices = @{&generate_slices_from_names($sa, \@slices, \@skip_slices, 'toplevel', $non_ref)};

#Doing this here for all the slices before we batch will be faster
my $rfsets = &get_regulatory_FeatureSets($analysis, \%cell_types);

#This will get rolled back again unnecessarily for each slice job, 
#but it should still be much quicker
#Well, it will be when rollback_FeatureSet supports multiple slices

if (! $local){ #BSUB!
  
  @tmp_args = @{&strip_param_flags(\@tmp_args, ('local', 'non_ref', 'include_mt'))};
  @tmp_args = @{strip_param_args(\@tmp_args, ('slices', 'skip_slices'))};

  (my $bsubhost = $host) =~ s/-/_/g;


  foreach my $slice(@slices){
	my $sr_name = $slice->seq_region_name;

	my $bsub = "bsub -q long -o $outdir/RegulatoryBuild_$sr_name.out -e $outdir/RegulatoryBuild_$sr_name.err ".
	  "-J 'RegulatoryBuild_$sr_name' -R 'select[my".$bsubhost."<80] rusage[my".$bsubhost."=10:duration=10]' -R 'select[mem>4000] rusage[mem=4000]' -M 4000000";
	
	my $cmd = "$ENV{EFG_SRC}/scripts/regulatory_build/build_regulatory_features.pl -local -slices $sr_name @tmp_args";
	  
	print "Submitting job for slice:\t".$sr_name."\n";
	system("$bsub $cmd") && die("Can't submit job array to farm");#Is this failure getting caught properly?

  }
  
  print "Output can be found here:\t$outdir\n";
  warn "We need a to put this in a Runnable so we can easily monitor which jobs have failed";

  exit;
}
elsif(scalar(@slices) > 1){
  die("Cannot run in local mode with more than one slice. Please omit -local or specify a single slice name.\n");
}


### Run slice job locally
my $slice = $slices[0];

# get core and fg seq_region_id for slice
my ($core_sr_id, $fg_sr_id);
$core_sr_id = $slice->get_seq_region_id;

# need to build cache first
$afa->build_seq_region_cache();
$fg_sr_id = $afa->get_seq_region_id_by_Slice($slice);
die("No eFG seq_region id defined for ".$slice->name." Almost certain there is no data".
	"available on this slice. You might want to run_update_DB_for_release.pl") 
  if ! defined $fg_sr_id;


### Dump AnnotatedFeatures
my $af_file = $outdir.'/annotated_features.'.$slice->seq_region_name.'.dat';

if ((! -f $af_file) ||
	(! $no_dump_annotated_features)){
  print "Dumping annotated features for slice:\t".$slice->name."\n";
  &dump_annotated_features();
} 
else {
  print "Using previously dumped annotated features:\t$af_file\n";
}

my $fh = open_file($af_file);


### Parse and build RegulatoryFeatures
if($debug_start && ! $debug_end){ #Set debug_end if not set
  $debug_end = $slice->end;
}

#Some variables which are used in subs
my (%afs, @rf, $ctype, $af_id, $sr_id, $sr_name, $start, $end, $strand, $score, $fset_id);
my (%feature_count, %seen_af, %removed_af);#Stats
my $rf_size = -1;

while (<$fh>) {
  # Read from file and process sequentially in sorted by start, end order. Each 
  # new feature is checked if it overlaps with the preceeding already seen 
  # features. If yes we just carry on with the next one. Otherwise
  next if (/^\#/);
  chomp;   
  ($af_id, $sr_id, $sr_name, $start, $end, $strand, $score, $fset_id) = split (/\s+/, $_);
  
  #Skip non debug regions
  #Warning, this may not build accurately!!!
  if ($debug_start) {
	next if ($start < $debug_start);
	last if ($start > $debug_end);
  }
    
  # Quick hack for 2nd/3rd version of reg. build. Need to disregard ChIPseq 
  # features below a certain threshold defined hardcoded in ChIPseq_cutoff hash.
  #next if (exists $ChIPseq_cutoff{$attrib_fsets{$fset_id}->name()} 
  #         && $score < $ChIPseq_cutoff{$attrib_fsets{$fset_id}->name()});
  
  print $_, "\t", $attrib_fsets{$fset_id}->name, "\n" if ($debug_start);
  my  $length = $end-$start+1;
  $ctype = $fset_cell_types{$fset_id};


#  if ($fset_id == 20){##
#
#   if($ctype ne 'CD4'){##
#
#	 die("CD4 ctype is $ctype");
#   }
#  }


 
  # some stats
  $feature_count{fsets}{$fset_id}++;
  $feature_count{ctypes}{$ctype}++; #Need to report these
	

  #if(($ctype eq 'IMR90') ||
  #	 ($ctype eq 'K562')){
  #	print "Seen $ctype\n";
  #  }

  if (exists $focus_fsets{$fset_id}) { ### Focus feature
	$feature_count{focus_ctypes}{$ctype}++;

	#Where is the overlap check for this?
	#Let's convert to multi-celltype first then fix
	
	if ($length > $focus_max_length) { # Too long, demote to attribute
	  print "Focus feature ($af_id) longer than $focus_max_length; ",
		"added to attribute feature list\n" if ($debug_start);

	  #Do we still want to build here even tho
	  #focus has exceeded the max length?

	  &add_feature();

	  #We are not adding this as an attribute feature
	  #Just removing it!!
	  
	  
	  $removed_af{$fset_id}{$af_id} = 1;
	  $removed_af{$ctype}++;
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
	
	#die("focus CD4 ctype is $ctype") if $fset_id = 149;

	&add_feature();
  } 
  else { ### Attribute feature
	
	#Only update if we have seen a focus feature of this ctype
	if( @rf && (exists ${$rf[$rf_size]{fsets}}{$ctype})){

	  if ($start <= $rf[$rf_size]{focus_end}){
		print "Adding overlapping $ctype feature(".$af_id.")\n" if ($debug_start);
	  
		# add annot. feature id to reg. feature
		$rf[$rf_size]{annotated}{$ctype}{$af_id} = undef;
		$rf[$rf_size]{fsets}{$ctype}{$fset_id}++;
		$seen_af{$fset_id}{$af_id} = 1;
		
		# update end of regulatory feature
		if ($end > $rf[$rf_size]{attribute_end}{$ctype}){
		  $rf[$rf_size]{attribute_end}{$ctype} = $end ;
		  print "Updating $ctype attr end to ".$end."\n" if $debug_start;
		}
	  } 
	  elsif( ($end <= $rf[$rf_size]{attribute_end}{$ctype}) &&
			 ($start <= ($rf[$rf_size]{focus_end} + $focus_extend)) ){		
		print "Adding fully contained $ctype feature(".$af_id.")\n" if ($debug_start);
		
		# add annot. feature id to reg. feature
		$rf[$rf_size]{annotated}{$ctype}{$af_id} = undef;
		$rf[$rf_size]{fsets}{$ctype}{$fset_id}++;
		$seen_af{$fset_id}{$af_id} = 1;
	  }
	}
	&add_feature();
  }
}


print "Skipping storage  RegulatoryFeatures" if ! $write_features;

&create_regulatory_features(\@rf);


if ($dump_regulatory_features) {
    
  #Change this to use @reg_feats from create_regulatory_features?

  my $outfile  = $outdir.'/regulatory_features_'.$slice->seq_region_name.'.dat';
  my $out = open_file($outfile, ">");
    
  map {
	printf $out "%d\t%s\t%d\t%d\t%d\t%d\t%s\t%s\n", 
	  $fg_sr_id, $slice->seq_region_name, 
        $_->{focus_start}, $_->{focus_end},
		  $_->{attribute_start}, $_->{attribute_end}, 
			$_->{binstring}, join(",", sort {$a<=>$b} keys %{$_->{annotated}});
  } @rf;
  
  printf "# Regulatory features written to ".$outfile."\n";
  
}



print "\n\n### RegulatoryBuild Report ###\n";

my $total_feature_count = 0;


#Do this in create_regulatory_seatures
#my (%rf_count);
#map { 
#  foreach my $k (keys %{$_->{fsets}}) {
#	#print join(" ", $k, $_->{fsets}->{$k}), "\n";
#	$rf_count{$k} += $_->{fsets}->{$k};
#  }
#} @rf;

my $die = 0;

foreach my $ctype(keys %ctype_fsets){
  next if $ctype eq 'core';

  printf "\n\n#$ctype FeatureSets %80s\t%8s\t%8s\t%8s\n",'Total', 'Included', 'Distinct', 'Removed';
  

  foreach my $set(@{$ctype_fsets{$ctype}}){
	my $set_feature_cnt = $feature_count{fsets}{$set->dbID} ||= 0;
	$total_feature_count += $set_feature_cnt;
	
  printf "#%80s(%d)\t%8d\t%8d\t%8d\t%8d\n", $set->name, $set->dbID, 
	$set_feature_cnt, #Total
	  $feature_count{regfeats}{$set->dbID} || 0, #Included
		scalar(keys %{$seen_af{$set->dbID}}) || 0, #Distinct
		  scalar(keys %{$removed_af{$set->dbID}}) || 0; #Removed
  }

   my $total = $feature_count{ctypes}->{$ctype} || 0;
  my $focus = $feature_count{focus_ctypes}->{$ctype} || 0;

  print "\n# Total $ctype features seen:\t$total\n";
  print "# Focus $ctype feature seen:\t$focus\n";

  if (($total < $focus) ||
	  $focus == 0){
	warn("Found inconsistenty between Total and Focus $ctype features\n");
	$die = 1;
  }

  #This number should be the Focus - those focus which have been removed
  my $regfs =  $feature_count{regfeats}{$ctype} || 0;
  print "# Total $ctype RegulatoryFeatures:\t$regfs\n";
  $removed_af{$ctype} ||= 0;
  my $expected_regfs = $regfs - $removed_af{$ctype};


  if( ($focus > $regfs) ||
	  ($expected_regfs != $regfs) ||
	  ($regfs == 0)){
	warn("Found inconsistency between Focus and true or expected($expected_regfs) RegulatoryFeature $ctype features\n");
	$die = 1;
  }
  

}


printf "# Total number of read features: %10d\n", $total_feature_count;
printf "\n# Total number of RegulatoryFeatures: %10d\n", scalar(@rf);

die("ERROR:\tInconsistencies were found between some of your focus sets and the resultant RegulatoryFeatures, see above.\n") if $die;

###############################################################################
# dump annotated features for given seq_region to file
sub dump_annotated_features () {
  #To do, remove this and just slurp directly into perl!

  my @fset_ids = keys %attrib_fsets;
  print "# Dumping ".scalar(@fset_ids)." AnnotatedFeature Sets for slice:\t".$slice->name."\n";

  my $sql = 'select annotated_feature_id, af.seq_region_id, sr.name,'.
	'seq_region_start, seq_region_end, seq_region_strand, score, feature_set_id '.
	  'from annotated_feature af, '.
		'(select distinct(seq_region_id), name from seq_region where seq_region_id='.$fg_sr_id.') sr '.
		  'where sr.seq_region_id=af.seq_region_id '.
			#"and schema_build='$schema_build' ".#This to stops the seq_region nr product
			#"and sr.name='".$slice->seq_region_name."' ".
			"and af.feature_set_id in (".join(',', @fset_ids).") order by seq_region_start, seq_region_end";
	
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
    
  print "# Execute: $command\n";

  # need to remove existing dump files, since we append to the file
  #system("rm -f $outdir/$dbname.annotated_feature_*.dat") &&
  #        die ("Can't remove files");
    
  #This is not catching mysql error due to the pipe.
  #Could do with counting features and validating vs number of lines

  system($command) && die ("Can't dump data to file in $outdir");
    
}

###############################################################################

#This should set up a new core RF
#We need to track annotation for all cell_type and then only build
#for those cell type for which we have we have seen core evidence
#We can do this by decoding the fset_ids

sub add_focus{
  print "\nNew $ctype focus\n" if ($debug_start);



  #if($ctype eq 'CD4'){
#	warn "CD4 $start $end";
#	#exit;
#  }

  push @rf, {
			 'focus_start' => $start,
			 'focus_end' => $end,
			 'attribute_start' => { $ctype => $start },
			 'attribute_end' => { $ctype => $end },
			 #focus features are migrated to annotated
			 #when feature is actually created
			 'focus' => { 
						 $ctype => { $af_id => undef },
						 core   => { $af_id => undef },
						},
			 'fsets' => {
						 core => {$fset_id => 1},
						 $ctype => {$fset_id => 1},
						}
			};

  #Do this here so we don't have to call $# multiple times.
  $rf_size +=1;  
  $seen_af{$fset_id}{$af_id} = 1;
  
  #print Dumper @rf;

}

sub update_focus{
  print "Updating focus(".join(', ', keys(%{$rf[$rf_size]{fsets}})).") with $ctype\n" if ($debug_start);

  #Set start values if we haven't seen this ctype before
  #update_attributes will be called immediately after this
  if(! exists  ${$rf[$rf_size]{fsets}}{$ctype}){
	print "Setting $ctype attr start/end to:\t$start - $end\n" if $debug_start;
	$rf[$rf_size]{attribute_start}{$ctype} = $start;
	$rf[$rf_size]{attribute_end}{$ctype}   = $end;
  }


  $rf[$rf_size]{focus_end} = $end if ($end > $rf[$rf_size]{focus_end});

 
  #Was always resetting attr end here, but we may have a longer feature already

  #Only set if previous cell_type attr was shorter
  if ($rf[$rf_size]{attribute_end}{$ctype} < $end) {
	print "Resetting $ctype attr end to:\t$end\n" if $debug_start;
	$rf[$rf_size]{attribute_end}{$ctype} = $end;
  }

  $rf[$rf_size]{fsets}{core}{$fset_id}++;
  $rf[$rf_size]{fsets}{$ctype}{$fset_id}++;
  $rf[$rf_size]{focus}{$ctype}{$af_id} = undef;
  $rf[$rf_size]{focus}{core}{$af_id} = undef;
  
  
  #Do we need to alter this to be celltype aware?
  $seen_af{$fset_id}{$af_id} = 1;

  #print Dumper @rf if ($debug_start);
  
}

sub update_attributes{
  print "Updating $ctype attributes\n" if ($debug_start);
  

  #This needs to be called only when we see a new ctype focus
  #Do we build on a ctype is there are encapsulated focus features
  #but non for a given ctype at this locus?
  
  # look upstream for all features overlapping focus feature


  for (my $i=0; $i<=$#{$afs{$ctype}}; $i++) {
	
	if ($afs{$ctype}->[$i]{end} >= $start) {
	  print "\tFound upstream overlapping attribute feature(".$afs{$ctype}->[$i]{af_id}.")\n" if $debug_start;
	  # FIRST update attribute start of the regulatory feature...
	  #Why is focus extend not being used here?

	  if ($afs{$ctype}->[$i]{start} < $rf[$rf_size]{attribute_start}{$ctype}) {
		#This always be the case the first feature
		#If there is one it will always be upstream or start at the same start as the new focus
		#And we process the rest in this block so no real need for this test?
		print "\tResetting $ctype attribute start to:\t".$afs{$ctype}->[$i]{start}."\n" if $debug_start;
		$rf[$rf_size]{attribute_start}{$ctype} = $afs{$ctype}->[$i]{start};
	  }
	  	  
	  # THEN remove features that are out of scope
	  splice(@{$afs{$ctype}}, 0, ($i));
		
	  #The current attr was not being added????
	  #attr end was not being set!!
	 
	  # NOW add the attribute features and set att end
	  map {
		if ($_->{end} >= $start ||
			$_->{start} >= $rf[$rf_size]{focus_start}-$focus_extend) {
		  
		  print "\tAdding attribute feature(".$_->{af_id}."):\t".$_->{start}.' - '.$_->{end}."\n" if $debug_start;

		  $rf[$rf_size]{annotated}{$ctype}{$_->{af_id}} = undef;
		  $rf[$rf_size]{fsets}{$ctype}{$_->{fset_id}}++;
		  
		  #Reset the attr end
		  if ($_->{end}   > $rf[$rf_size]{attribute_end}{$ctype}) {
			print "\tResetting $ctype attribute end to:\t\t".$_->{end}."\n" if $debug_start;
			$rf[$rf_size]{attribute_end}{$ctype} = $_->{end};
		  }
		  # could also test for demoted focus features here!
		  
		  $seen_af{$fset_id}{$af_id} = 1;
		}
	  } @{$afs{$ctype}};
	  
	  last;
	}
  }
}

sub add_feature{
  
  #print "Adding $ctype feature($af_id) to attribute feature_list\n" if $debug_start;


  push(@{$afs{$ctype}}, 
	   {
		af_id => $af_id,
		start => $start,
		end => $end,
		strand => $strand,
		score => $score,
		fset_id => $fset_id
	   });
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
  
  #print Dumper $rf;
  #print Dumper $rfset->feature_type;
  #print $rfset->feature_type;

  print "\n\nCreating ".scalar(@$rfs)." RegulatoryFeatures\n" if $debug_start;

  foreach my $rf (@$rfs) {
	
	if ($slice->start != 1 || $slice->strand != 1) {
	  warn("**** SLICE doesn't start at position 1; resetting start of slice to 1 ****");
	  $slice->{'start'} = 1;
	  #This only works as we are working with seq_region_start/end values in the dumps
	  #We should really fetch the slices with the include_duplcates flag i.e. 1 as start for PAR regions
	}

	#Now build each ctype RegFeat
	foreach my $ct (keys %{$rf->{fsets}}) {

	  #Only count this is it is a true regf and not a scaffold projection
	  $feature_count{regfeats}{$ct}++;
	  
	  map $feature_count{regfeats}{$_} += $rf->{fsets}{$ct}{$_}, keys %{$rf->{fsets}{$ct}};#$_ is fset_id

	  #Build attr cache dependant on presence of attrs
	  my $attr_cache;

	  if(exists $rf->{annotated}->{$ct}){
		$attr_cache = {
					   %{$rf->{focus}->{$ct}},
					   %{$rf->{annotated}->{$ct}},
					  };
	  }
	  else{
		$attr_cache = {
					   %{$rf->{focus}->{$ct}},
					  };
	  }

	  print "\nBuilding $ct RegulatoryFeature with attribute cache:\n".Data::Dumper::Dumper($attr_cache)."\n" if $debug_start;

	  my $reg_feat = Bio::EnsEMBL::Funcgen::RegulatoryFeature->new
		(
		 -slice            => $slice,
		 -start            => $rf->{focus_start},
		 -end              => $rf->{focus_end},
		 -strand           => 0,
		 -display_label    => &build_binstring($rf, $ct),
		 -feature_set      => $rfsets->{$ct},
		 -feature_type     => $rfsets->{$ct}->feature_type,
		 -_attribute_cache => {'annotated_feature' => $attr_cache },
		);
	  
	  if ($write_features) {	#Store
		  
		($reg_feat) = @{$rfa->store($reg_feat)};
		  
		#Sanity check bound here to make sure we have built correctly
		#Can only test bounds if adaptors are set

		if($ct ne 'core'){#Only have bound for non-core features

		  if (($reg_feat->bound_start != $rf->{attribute_start}->{$ct}) ||
			  ($reg_feat->bound_end != $rf->{attribute_end}->{$ct})) {
			
			warn $reg_feat->bound_start.' != '.$rf->{attribute_start}->{$ct};
			warn $reg_feat->bound_end .' != '.$rf->{attribute_end}->{$ct};
			die('Calculated attribute_start/end values do not match bound_start/end values');
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
  warn "Need to implement Helper::define_and_validate_sets here";
  #Can we not use Helper::define_and_validate_sets for this?
  #Just hack for now!
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
	

  foreach my $ctype (keys %{$ctypes}) {	
	my ($desc, $fset_name);
	
	if($ctype eq 'core'){
	  $fset_name = 'RegulatoryFeatures';
	  $desc = 'Generic RegulatoryFeature focus regions';
	}
	else{
	  $fset_name = "RegulatoryFeatures:$ctype";
	  $desc = "$ctype specific RegulatoryFeatures";
	  
	}

	my $rollback = ($clobber) ? 'supporting_sets' : undef;

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
												 \@slices);

	$rf_sets{$ctype} = $dset->product_FeatureSet;

  }
    
  $helper->log("Got RegulatoryFeature sets for CellTypes:\t".join(', ', keys %rf_sets)); 
  return \%rf_sets;
}

1;
