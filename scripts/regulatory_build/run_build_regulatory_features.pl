#!/software/bin/perl
###!/usr/bin/perl

=head1 NAME

run_build_regulatory_features.pl -- wrapper script to run the
build_regulatory_features.pl to create "Ensembl Regulatory Build"

=head1 SYNOPSIS

run_build_regulatory_features.pl -host host -user user -pass password 


Options:

  Mandatory
    -host|h          Host for eFG DB
    -user|u          User for eFG DB
    -pass|p          Password for eFG DB
    -dbname|d        Name of eFG DB
    -data_version|v  Version of data in eFG DB (e.g. 51_36m)
    -outdir|o        Name of outputut directory

    -cdbhost
    -cdbport
    -cdbuser
    -cdbpass

  Optional
    -jobs
    -help|?
    -man|m

=head1 DESCRIPTION

Wrapper script to handles build job submissions, set generation, rollback and other 
high level functions.

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
#1 Integrate with build_regulator_features.pl and just have two modes, one to do the set validation, dumping
#and submitting, and the other to actually run the jobs??? Are these now sufficiently different ot keep separate?

#2 Update to take standard pipeline.env params i.e. remove dataversion and define dnadb explicitly.

#3 Remove dump_annotated_feature func from build_regulatory_features and submit parallelised jobs
#to farm using export_features script/Exporter.pm

#4 Store data_set!?

#5 Check running jobs

#6 DONE Add archive function!

use strict;
use warnings;
use Data::Dumper;

use Bio::EnsEMBL::Utils::Exception qw(verbose warning info);
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::Helper;

use Getopt::Long;

my ($dbhost,$dbport,$dbuser,$dbpass,$dbname,$species, $fset, $outdir, $seq_region_name,
	$cdbhost,$cdbport,$cdbuser,$cdbpass,$cdbname, $non_ref, $include_mt, $version,
	@focus_names, @attr_names, $cdb, $db, $jobs,$help,$man, $dump_afs);

$|=1;

# default core database 
$cdbhost ||= 'ens-staging';
$cdbport ||= 3306;
$cdbuser ||= 'ensro';


GetOptions (
			#Use pipeline env db param names
            "dbhost|h=s"       => \$dbhost,
            "dbport=s"         => \$dbport,
            "dbuser|u=s"       => \$dbuser,
            "dbpass|p=s"       => \$dbpass,
            "dbname|d=s"     => \$dbname,
            "species=s"      => \$species,
			"dnadb_host=s"      => \$cdbhost,
            "dnadb_port=s"      => \$cdbport,
            "dnadb_user=s"      => \$cdbuser,
            "dnadb_pass=s"      => \$cdbpass,
			'dnadb_name=s'   => \$cdbname,
			'seq_region_name=s'   => \$seq_region_name,
			'focus_sets=s{,}'     => \@focus_names, 
			'attribute_sets=s{,}' => \@attr_names,
			'version=s'           => \$version,
			'outdir=s'            => \$outdir,
			'non_ref'             => \$non_ref,
			'include_mt'          => \$include_mt,
			'dump_annotated_features' => \$dump_afs,
			
			#To add opts
			#slice
			#skip_slices


            #"jobs=s"         => \$jobs,
            "help|?"         => \$help,
            "man|m"          => \$man,
            );

### check options ###

die("Must specify mandatory database hostname (-dbhost).\n") if ! defined $dbhost;
die("Must specify mandatory database username. (-dbuser)\n") if ! defined $dbuser;
die("Must specify mandatory database password (-dbpass).\n") if ! defined $dbpass;
die("Must specify mandatory database name (-dbname).\n") if ! defined $dbname;

#Need to validate whe have sets pass

#Allow comma separated
@focus_names = split(/,/,join(',',@focus_names));
@attr_names = split(/,/,join(',',@attr_names));


if($cdbname){
  $cdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (
     -host => $cdbhost,
     -port => $cdbport,
     -user => $cdbuser,
     -dbname => $cdbname,
     -species => $species,
	 -group => 'core',
     );
}


#Can we pass just dnadb_host here, to avoid having to sepcify name in run script?
$db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
    (
     -host   => $dbhost,
     -user   => $dbuser,
     -dbname => $dbname,
     -species => $species,
     -pass   => $dbpass,
     -port   => $dbport,
     -dnadb  => $cdb,
	 -group => 'funcgen',
     );



### CHECK ARCHIVED VERSIONS
#so we don't clobber it before we have a chance to stable_id map
#Or should we remove this facility of host more than one reg build and map between DBs?
#Move this to Healthchecker for reuse in update_DB_for_release?

if(! $version){
  die('To properly archive a previous regulatory build, you must provide the version number of the new build');
}

my $old_version = $version -1;
my $mc = $db->get_MetaContainer;
my $some_old_not_archived = 0;
my $some_old_archived = 0;
my $fsa = $db->get_FeatureSetAdaptor();
my @meta_keys = ('regbuild.feature_set_ids', 'regbuild.initial_release_date', 'regbuild.feature_type_ids', 'regbuild.last_annotation_update');

#Now check meta entries for old build
my ($current_version) = @{$mc->list_value_by_key('regbuild.version')};

if((! defined $current_version) && 
   ($version != 1)){
  die('Could not find regbuild.version meta entry, please correct manually');
}

#This is wrong as we may have already updated and just want to clobber
if(($current_version != $old_version) &&
   ($current_version != $version)){
  die("The regbuild version specified($version) is not the next or current version according to meta regbuild.version($current_version). Please correct manually\n");
}


for my $mkey(@meta_keys){
  #Check both the current and the old versions just in csae we forgot to update the version when running

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
else{ #test old and new feature set here for sanity
  #So we aren't depending solely on the meta keys
  my $fset = $fsa->fetch_by_name("RegulatoryFeatures_v${version}");

  if($fset){
	die("It appears that the new RegulatoryBuild feature set(v${version}) has already been archived, but no meta entries are present. Please correct manually.\n")
  }

  if($version != 1){
	$fset = $fsa->fetch_by_name("RegulatoryFeatures_v${old_version}");
  }
 
  if($some_old_not_archived){#all old meta not archived
	
	if($fset){
	  die("It appears that the old RegulatoryBuild feature set(v${old_version}) has already been archived, but no meta entries are present. Please correct manually.\n")
	}

	if($version != 1){
	  print "Achiving old RegulatoryFeatures set\n";
	 
	
	  #We shouldn't need to track this process with a status
	  #As we have the exhaustive tests above.
	  
	  #Rollback all other older feature_sets?

	  #rename feature_set
	  my $sql = "update feature_set set name='RegulatoryFeatures_v${old_version}' where name='RegulatoryFeatures'";
	  $db->dbc->db_handle->do($sql);

	  #validate update and revoke states
	  my $fset = $fsa->fetch_by_name("RegulatoryFeatures_v${old_version}");

	  if(! $fset){
		die("Failed to create archive FeatureSet:\tRegulatoryFeatures_v${old_version}\n");
	  }
	  
	  $fset->adaptor->revoke_status('MART_DISPLAYABLE', $fset);
	  $fset->adaptor->revoke_status('DISPLAYABLE', $fset);
	  #Keep DAS_DISPLAYABLE status?

	  #rename meta keys
	  foreach my $mkey(@meta_keys){
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
	  die("It appears that all the old RegulatoryBuild meta entries have been archived, but no archived feature set(v${old_version}) is present. Please correct manually.\n")
	}

	print "Old RegulatoryFeature set has previously been archived\n";

  }
}



# determine the number of toplevel slices we can analyze 
# (excluding mitochondria sequences (MT)!!!)
# this bit needs to be moved t the pipeline helper that it also can be called 
# by the actual run script

my @slices = ();
my $sa = $db->dnadb->get_SliceAdaptor();



#Change this back it was a problem with the indices

if($seq_region_name){
  

  my $slice = $sa->fetch_by_name($seq_region_name);

  if(! $slice){
	die("You have specified a slice which is not in the db:\t$seq_region_name\n");
  }
  
  @slices = ($slice);

  print "Running with slice:\t".$slice->name."\n";
}
else{

  foreach my $s (sort {$a->name cmp $b->name} @{$sa->fetch_all('toplevel', undef, $non_ref)}){
	
	next if ($s->seq_region_name eq 'MT' && (! $include_mt));
	push @slices, $s;
	
  }
  
  $non_ref = "non_ref" if $non_ref;
  print 'Running with '.scalar(@slices)." toplevel $non_ref slices\n"
}


#Don't need to batch this, there will probably be less than 30!
#$jobs = '1-'.scalar(@slices) unless (defined $jobs);




#Much of this may be repeated in the build script for validation?

# parse both focus and attribute sets and check that they exist
my (%focus_fsets, %attrib_fsets);


map { $fset = $fsa->fetch_by_name($_);
      die("Focus set $_ does not exist in the DB") 
          if (! defined $fset); 
      $focus_fsets{$fset->dbID} = $fset; 
  } @focus_names;
#print Dumper %focus_fsets;

map { 
    $fset = $fsa->fetch_by_name($_);
    die("Attribute set $_ does not exist in the DB") 
        if (! defined $fset); 
    $attrib_fsets{$fset->dbID()} = $fset; 
} @attr_names;
#print Dumper %attrib_fsets;

# make sure that attribute sets also contain focus sets (Do we really need this?)
map { $attrib_fsets{$_} = $focus_fsets{$_} } keys %focus_fsets;


my (@fset_ids, @ftype_ids);

foreach $fset(sort {$a->name cmp $b->name} values %attrib_fsets) {
  push @fset_ids, $fset->dbID;
  push @ftype_ids, $fset->feature_type->dbID;
}


print "Fset string is:\t".join(',', @fset_ids)."\n";
print "Ftype string is:\t".join(',', @ftype_ids)."\n";

my ($sql, $meta_value, $reg_string, $cmd);


### build and import regbuild strings by feature_set_id and feature_type_id

my $sth = $db->dbc->prepare("select meta_value from meta where meta_key='regbuild.feature_set_ids'");
$sth->execute();
($meta_value) = map "@$_", @{$sth->fetchall_arrayref};



$reg_string = join(',', map {$_->dbID} sort {$a->name cmp $b->name} values %attrib_fsets);

if (! $meta_value) {


    $sql = "insert into meta (meta_key, meta_value) values ('regbuild.feature_set_ids', '$reg_string')";

    eval {
        $db->dbc->do($sql);
    };
    die("Couldn't store regbuild.feature_set_ids in meta table.\n$@") if ($@);

} else {

  if($meta_value ne $reg_string){
	die("regbuild.feature_set_ids already exists and does not match\nOld\t$meta_value\nNew\t$reg_string\nPlease archive previous RegulatoryBuild.\n");
  }else{  
	warn "Found matching regbuild.feature_set_ids meta entry, do you need to archive a previous Regulatorybuild?\n";
	#Skip this warn with clobber?
	#Add this warn to summary log?

  }
}


$sth = $db->dbc->prepare("select meta_value from meta where meta_key='regbuild.feature_type_ids'");
$sth->execute();
($meta_value) =  map "@$_", @{$sth->fetchall_arrayref};
$reg_string = join(',', map {$_->feature_type->dbID} sort {$a->name cmp $b->name} values %attrib_fsets);

if (! $reg_string) {

    $sql = "insert into meta (meta_key, meta_value) values ('regbuild.feature_type_ids', '$reg_string')";

    eval {
        $db->dbc->do($sql);
    };
    die("Couldn't store regbuild.feature_type_ids in meta table.") if ($@);

} else {
 if($meta_value ne $reg_string){
	die("regbuild.feature_type_ids already exists and do not match\nOld\t$meta_value\nNew\t$reg_string\nPlease archive previous RegulatoryBuild.\n");
  }else{  
	warn "Found matching regbuild.feature_type_ids meta entry, do you need to archive a previous Regulatorybuild?\n";

	#Skip this warn with clobber?
	#Add this warn to summary log?

  }    
}

#Enough time to realise we might want to archive
sleep(5);



# submit jobs to the farm

(my $bsubhost = $dbhost) =~ s/-/_/g;

$dump_afs = ' -dump_annotated_features ' if $dump_afs;

foreach my $slice(@slices){
  my $sr_name = $slice->seq_region_name;
  

  #my $bsub = "bsub -oo '$outdir/RegulatoryBuild_%I.log' -J 'RegulatoryBuild_[".$jobs."]' ".
  #  "-R 'select[my".$bsubhost."<80] rusage[my".$bsubhost."=10:duration=10]'";

 my $bsub = "bsub -o $outdir/RegulatoryBuild_$sr_name.log -e $outdir/RegulatoryBuild_$sr_name.err -J 'RegulatoryBuild_$sr_name' ".
   "-R 'select[my".$bsubhost."<80] rusage[my".$bsubhost."=10:duration=10]'";
 

  $cmd = "$ENV{EFG_SRC}/scripts/regulatory_build/build_regulatory_features.pl"
	.' -seq_region_name '.$slice->name
    .' -host '.$dbhost
	  .' -port '.$dbport
		.' -user '.$dbuser
		  .' -pass '.$dbpass
			.' -dbname '.$dbname
			  #.' -species '.$species
				.' -dnadb_host '.$cdbhost
				  .' -dnadb_name '.$cdbname
					.' -focus '.join(',', @focus_names)
					  .' -attrib '.join(',', @attr_names)
						.' -outdir '.$outdir
						  .' -dump_regulatory_features '
							.' -clobber '
							  .' -stats '
								." -write_features $dump_afs";

  #print $bsub, "    ", $cmd, "\n";
  
  #system($cmd);
  #exit;

  system("$bsub $cmd") && die("Can't submit job array to farm");
}

