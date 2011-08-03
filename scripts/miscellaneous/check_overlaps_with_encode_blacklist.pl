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

check_overlaps_with_blacklist.pl - For all the human annotated features, checks how they overlap with a set of locations defined as being problematic by the ENCODE project (so called blacklist). Produces log files of the overlap analysis.

=head1 SYNOPSIS

check_overlaps_with_blacklist.pl -efgdb_host dbhost -efgdb_port dbport -efgdb_user dbuser -efgdb_name dbname\
    -coredb_host dnadbhost -coredb_port dnadbport -coredb_user dnadbuser -coredb_name dnadbname\
    -remove -output overlaps.txt -feature_set fset1 fset2 ...

Options:

  Mandatory:
    -feature_sets        Space-separated list of names of feature sets for which to check the overlap

  Optional:
    -efgdb_host          Host for eFG DB [$DB_HOST]
    -efgdb_port          Port for eFG DB [$DB_PORT]
    -efgdb_user          User for eFG DB [$DB_USER]
    -efgdb_name          Name of eFG DB [$DB_NAME]
    -efgdb_pass          Password fof the eFG DB
    -coredb_host         Host for Core DB [$DNADB_HOST]
    -coredb_port         Port for Core DB [$DNADB_PORT]
    -coredb_user         User for Core DB [$DNADB_USER]
    -coredb_name         Name of Core DB [$DNADB_NAME]
    -remove              When specified it removes the overlapping features from the database 
    -output              File to store results of analysis  [overlaps.txt]
    -all                 When specified all sets are considered

=head1 DESCRIPTION

Downloads datasets in the data tracking database. Updates their download status.

=head1 OPTIONS

=over

=item B<help>

Gives this help menu

=item B<-all>

When specified, all sets are considered

=item B<-remove>

When specified, it removes features overlapping with encode blacklisted regions (requires database user with write privileges)

=item B<-output>

Name of output file to create (defaults to output.txt)

=item B<-coredb_host>

Host where the core database is (defaults to $DNADB_HOST)

=item B<-coredb_user>

User of the core database (defaults to $DNADB_USER)

=item B<-coredb_pass>

Password for the core database user (defaults to $DNADB_PASS)

=item B<-coredb_port>

Port of the host where the core database is (defaults to $DNADB_PORT)

=item B<-coredb_name>

Name of the Core database (defaults to $DNADB_NAME)

=item B<-efgdb_host>

Host where the EFG database is (defaults to $DB_HOST)

=item B<-efgdb_user>

User of the EFG database (defaults to $DB_USER)

=item B<-efgdb_pass>

Password for the EFG database user (defaults to $DB_PASS)

=item B<-efgdb_port>

Port of the host where the EFG database is (defaults to $DB_PORT)

=item B<-efgd_bname>

Name of the EFG database (defaults to $DB_NAME)

=back

=cut


#To do
# 1 Set BLACKLIST_FILTERED status?
# 2 add -slices -skip_slices support
# 3 Probably want to implement this in a module so it can be re-used in the pipeline
# 4 refine log output

use warnings;
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);

#Variables from the EFG and pipeline environments
my $efgdb_host = $ENV{DB_HOST};
my $efgdb_port = $ENV{DB_PORT};
my $efgdb_user = $ENV{DB_USER};
my $efgdb_name = $ENV{DB_NAME};
my $efgdb_pass;

my $coredb_host = $ENV{DNADB_HOST};
my $coredb_port = $ENV{DNADB_PORT};
my $coredb_user = $ENV{DNADB_USER};
my $coredb_name = $ENV{DNADB_NAME};

my $species = 'homo_sapiens';
my $output = "overlaps.txt";
my $all;
my $remove;
my $help;
my @feature_sets;

GetOptions ("output=s"      => \$output,
	    "efgdb_port=s"  => \$efgdb_port,
            "efgdb_host=s"  => \$efgdb_host,
            "efgdb_user=s"  => \$efgdb_user,
            "efgdb_name=s"  => \$efgdb_name,
            "efgdb_pass=s" => \$efgdb_pass,
            "coredb_port=s" => \$coredb_port,
            "coredb_host=s" => \$coredb_host,
            "coredb_user=s" => \$coredb_user,
            "coredb_name=s" => \$coredb_name,
			"all"           => \$all,
			#"class"         => \$feature_class,
			#Can't implement other classses here as they have xrefs etc
			#would have to use Helper to rollback feature methods
            "remove"        => \$remove,
            "feature_sets=s{,}"  => \@feature_sets,
	    "help|h"              => \$help,
	   )  or pod2usage( -exitval => 1 ); #Catch unknown opts

pod2usage(1) if ($help || !$efgdb_name);

#Get db adaptors: 

my $coredba = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (
     -host => $coredb_host,
     -port => $coredb_port,
     -user => $coredb_user,
     -dbname => $coredb_name,
     -species => $species,
     -group   => 'core',
     );

if(!$coredba){ warn "Could not connect to core database..."; exit 1; }

my $efgdba = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
    (
     -host   => $efgdb_host,
     -user   => $efgdb_user,
     -pass   => $efgdb_pass,
     -dbname => $efgdb_name,
     -species => $species,
     -port   => $efgdb_port,
     -dnadb  => $coredba,
     -group  => 'funcgen',
     );

if(!$efgdba){ warn "Could not connect to EFG database..."; exit 1; }

my $fsa = $efgdba->get_FeatureSetAdaptor();
my $afa = $efgdba->get_AnnotatedFeatureAdaptor();

#my $cta = $efgdba->get_CellTypeAdaptor();
#my $K562 = $cta->fetch_by_name("K562");
#This doesn't seem to be working very well... TODO: check why...
#my @fsets = $fsa->fetch_all_by_CellType($K562);


#Currently only handles annotated FeatureSets!

if((scalar(@feature_sets)>0) && $all){ warn " Only specified feature sets will be considered. -all is being ignored."; }

my @fsets;

if(! $all){
  foreach my $feature_set (@feature_sets){ 
	my $fset = $fsa->fetch_by_name($feature_set);
	
	if($fset) { 
	  push(@fsets, $fsa->fetch_by_name($feature_set)); 
	} 
	else { warn "Could not find Feature Set $feature_set"; }
  }
}
else{
  
  #warn "Features for external or regulatory sets overlapping blacklist are not removed, only counts are reported!";
  #@fsets = @{$fsa->fetch_all()};
  #removed this as it was throwing errors when trying to fetch via the specific feature adaptor
  #would have to specify correct feature adaptors

  @fsets = @{$fsa->fetch_all_by_type('annotated')};
  warn "All Annotated FeatureSets(".(scalar(@fsets)).") are being used... this may take a while!";
  
}


if(scalar(@fsets) ==0){
  die("No Feature Sets found. Use -feature_sets or -all");
}

my %blacklist;
open(FO,">".$output) || die("Cannot open:\t$output");
my @count_keys = ('count', 'total_feature_length', 'cumulative_region_coverage');
print FO "FeatureSet\t".join("\t", @count_keys)."\n";


#add global counts here?
my %global_counts;
my $slice_a = $coredba->get_SliceAdaptor();
#Include duplicate regions...

foreach my $chromosome_slice (@{ $slice_a->fetch_all('toplevel', undef, 0, 1) }){
  

  #my @set_features = @{$fset->get_Features_by_Slice($chromosome_slice)};
  #Do we really want to do this set wise here or use af_a->fetch_all_by_Slice_FeatureSets

  
  my @excluded_regions = map $_->feature_Slice, 
	@{$chromosome_slice->get_all_MiscFeatures('encode_excluded')}; 
  
  if ($chromosome_slice->seq_region_name eq 'Y') {
	print "Adding the Y PAR regions to the blacklist\n";
	#We need to add them at the beginning so that these large regions will be tested first... 
	unshift @excluded_regions, $coredba->get_SliceAdaptor()->fetch_by_region('chromosome','Y',10001,2649521);
	unshift @excluded_regions, $coredba->get_SliceAdaptor()->fetch_by_region('chromosome','Y',59034050);
  }


  foreach my $region (@excluded_regions) {
	my @ids_to_remove;
	my %counts;

	#could probably do with some global counts too.


	print "Fetching features from black list region:\t\t\t\t".$region->name."\n";
	#my @set_features = @{$fset->get_Features_by_Slice($regions)};
	my @set_features = @{$afa->fetch_all_by_Slice_FeatureSets($region, \@fsets)};
	
	
	#Could use has_par here to set sr_id_test flag
	#do we even need the feature loop?
	#Yes, if we want to maintain the overlap stats and %blacklist
	my $has_par = ($region->seq_region_name eq 'Y') ? 1 : 0;
	
	#Need to overhaul stats/logs here
	#but should work

	foreach my $feature (@set_features) {
		
	  #Firstly avoid PAR regions which have been projected via the slice fetch code
	  if($has_par){
		my $query_sr_id = $slice_a->get_seq_region_id($region);
		my $sql = 'select seq_region_id from annotated_feature where annotated_feature_id='.$feature->dbID;
		my ($db_sr_id) = $efgdba->dbc->db_handle->selectrow_array($sql);
		next if $query_sr_id != $db_sr_id;
	  }

	  $global_counts{$feature->feature_set->name}{count} ++;
	  $counts{$feature->feature_set->name}{count} ++;
	  $counts{$feature->feature_set->name}{total_feature_length} += $feature->feature_Slice()->length();


	  
	  #This assumes that encode_excluded regions are non-overlapping with each other
	  #There is an overlap
	  # This has the disadvantage that overlap counts may be more than the features... 
	  #$set_overlap_count++;  
	  push @ids_to_remove, $feature->dbID();

	  my $start = $region->start();
	  my $end = $region->end();
	  
	  if ($start<$feature->start()) { 
		$start = $feature->start();  
	  }
		
	  if ($end>$feature->end()) { 
		$end = $feature->end();  
	  }
		
	  $counts{$feature->feature_set->name}{cumulative_region_coverage} += ($end-$start+1);

		
	  #The same feature set can overlap several times with the same blacklist region...
	  #Can probably update this or integrate with other log?

	  #What is this for?

	  $blacklist{$region->name}{$feature->feature_set->name} = 1;
	  
	}
  

	foreach my $fset_name(keys %counts){

	  print FO "${fset_name}\t".join("\t", (map $counts{$fset_name}{$_}, @count_keys))."\n";
	}
    
	&remove_features($region, \@ids_to_remove);
  }
}


print FO "\n\nTotal counts\n";

foreach my $fset_name(keys %global_counts){

  print FO "${fset_name}\t".$global_counts{$fset_name}."\n";
}


close FO;


open(FOB,">".$output.".blacklist");
foreach my $region (sort keys %blacklist){
  print FOB $region;
  foreach my $set (keys %{$blacklist{$region}}){ print FOB "\t".$set; }
  print FOB "\n";
}
close FOB;


sub remove_features{
  my ($region, $db_ids, $counts) = @_;

  
  
  if(scalar(@$db_ids) == 0){
	print "No AnnotatedFeatures for blacklist region:\t".$region->name."\n";
  }
  elsif($remove){
	#TODO Batch remove features!
		print "Deleting ".scalar(@$db_ids).
		  " AnnotatedFeatures for blacklist region:\t".$region->name."\n";

	my $sql = "DELETE FROM annotated_feature WHERE annotated_feature_id in (".join(', ', @$db_ids).')';
	
	eval{ $efgdba->dbc->do($sql);	};
	
	if($@){ 
	  die("Could not delete ".scalar(@$db_ids).
		  " AnnotatedFeatures for blacklist region:\t".$region->name."\n$@");
	}
	
	# foreach my $fID (keys %ids_to_remove){
	#   eval{ 
	#     my $sql = "DELETE FROM annotated_feature WHERE annotated_feature_id=".$fID.";";
	#     $efgdba->dbc->do($sql);		 
	#   };
	#   if($@) { warn $@->getErrorMessage(); }
  }
  else{
	print "Skipping deletion of ".scalar(@$db_ids).
	  " AnnotatedFeatures for blacklist region:\t".$region->name."\n";
  }
}




1;


