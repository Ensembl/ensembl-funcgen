#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 NAME

check_overlaps_with_blacklist.pl - For all the human annotated features, checks how they overlap with a set of locations defined as being problematic by the ENCODE project (so called blacklist). Produces log files of the overlap analysis.

=head1 SYNOPSIS

check_overlaps_with_blacklist.pl -dbhost dbhost -dbport dbport -dbuser dbuser -dbname dbname\
    -dnadb_host dnadbhost -dnadb_port dnadbport -dnadb_user dnadbuser -dnadb_name dnadbname\
    -remove -output overlaps.txt -feature_set fset1 fset2 ...

Options:

  Mandatory:
    -feature_sets        Space-separated list of names of feature sets for which to check the overlap

  Optional:
    -dbhost          Host for eFG DB [$DB_HOST]
    -dbport          Port for eFG DB [$DB_PORT]
    -dbuser          User for eFG DB [$DB_USER]
    -dbname          Name of eFG DB [$DB_NAME]
    -dbpass          Password fof the eFG DB
    -dnadb_host      Host for Core DB [$DNADB_HOST]
    -dnadb_port      Port for Core DB [$DNADB_PORT]
    -dnadb_user      User for Core DB [$DNADB_USER]
    -dnadb_name      Name of Core DB [$DNADB_NAME]
    -remove          When specified it removes the overlapping features from the database 
    -output          File to store results of analysis  [overlaps.txt]
    -all             When specified all sets are considered

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

=item B<-dnadb_host>

Host where the core database is (defaults to $DNADB_HOST)

=item B<-dnadb_user>

User of the core database (defaults to $DNADB_USER)

=item B<-dnadb_pass>

Password for the core database user (defaults to $DNADB_PASS)

=item B<-dnadb_port>

Port of the host where the core database is (defaults to $DNADB_PORT)

=item B<-dnadb_name>

Name of the Core database (defaults to $DNADB_NAME)

=item B<-dbhost>

Host where the EFG database is (defaults to $DB_HOST)

=item B<-dbuser>

User of the EFG database (defaults to $DB_USER)

=item B<-dbpass>

Password for the EFG database user (defaults to $DB_PASS)

=item B<-dbport>

Port of the host where the EFG database is (defaults to $DB_PORT)

=item B<-dbname>

Name of the EFG database (defaults to $DB_NAME)

=back

=cut


#To do
# 1 Set BLACKLIST_FILTERED status?
# 2 add -slices -skip_slices support
# 3 Move to a module so it can be re-used in the pipeline and via stand alone script
# 4 refine log output
# 5 add RF and MF support? Shouldn't need this is we always filter after peak calling.

use warnings;
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);

#Variables from the EFG and pipeline environments
my ($efgdb_host, $efgdb_port, $efgdb_user, $efgdb_name, $efgdb_pass);
my ($coredb_host, $coredb_port, $coredb_user, $coredb_name);
my $output = "overlaps.txt";
my ($chr, $all, $remove, $help, @feature_sets);

GetOptions
  (
   "output=s"     => \$output,
   "dbport=s"     => \$efgdb_port,
   "dbhost=s"     => \$efgdb_host,
   "dbuser=s"     => \$efgdb_user,
   "dbname=s"     => \$efgdb_name,
   "dbpass=s"     => \$efgdb_pass,
   "dnadb_port=s" => \$coredb_port,
   "dnadb_host=s" => \$coredb_host,
   "dnadb_user=s" => \$coredb_user,
   "dnadb_name=s" => \$coredb_name,
   "all"          => \$all,
   #"class"       => \$feature_class,
   #Can't implement other classses here as they have xrefs etc
   #would have to use Helper to rollback feature methods
   "remove"       => \$remove,
   'chr=s'        => \$chr,
   "feature_sets=s{,}"  => \@feature_sets,
   "help|h"              => \$help,
  )  
  or pod2usage( -exitval => 1 ); #Catch unknown opts

pod2usage(1) if ($help || !$efgdb_name);

#Get db adaptor 
my $efgdba = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
    (
     -host   => $efgdb_host,
     -user   => $efgdb_user,
     -pass   => $efgdb_pass,
     -dbname => $efgdb_name,
     -port   => $efgdb_port,
     -dnadb_host => $coredb_host,
     -dnadb_port => $coredb_port,
     -dnadb_user => $coredb_user,
     -dnadb_name => $coredb_name,
     );

if($efgdba->species ne 'homo_sapiens'){
  die('This script currently only works for homo_sapiens');
}

if(! defined $efgdba){ 
  die("Could not connect to the $efgdb_name database");
}

my $coredba = $efgdba->dnadb;

#Check connections
$efgdba->dbc->db_handle;
$coredba->dbc->db_handle;


my $fsa     = $efgdba->get_FeatureSetAdaptor;
my $afa     = $efgdba->get_AnnotatedFeatureAdaptor;


if((scalar(@feature_sets)>0) && $all){
  die("You have specified mutually exclusive parameters, please -all or -feature_sets"); 
}

my ($fset, @fsets);

if(! defined $all){
  
  foreach my $feature_set (@feature_sets){ 
    $fset = $fsa->fetch_by_name($feature_set);
    
    if($fset) { 
      push @fsets, $fsa->fetch_by_name($feature_set); 
    } 
    else { 
      die("Could not find Feature Set:\t$feature_set"); 
    }
  }
}
else{
  
  #warn "Features for external or regulatory sets overlapping blacklist are not removed, only counts are reported!";
  #@fsets = @{$fsa->fetch_all()};
  #removed this as it was throwing errors when trying to fetch via the specific feature adaptor
  #would have to specify correct feature adaptors
  @fsets = @{$fsa->fetch_all_by_type('annotated')};
}

if(scalar(@fsets) ==0){
  die("No Feature Sets found");
}

my (%blacklist, %global_counts);
#These counts do not contain feature sets which have no
#filtered features

my $slice_a = $coredba->get_SliceAdaptor;

open(FO,">".$output) || die("Cannot open:\t$output");
#Force order with array here for printing
my @count_keys = ('count', 'total_feature_length', 'cumulative_region_coverage');
print FO "FeatureSet\t".join("\t", @count_keys)."\n";


### Fetch all toplevel slices
#Include duplicate regions (we want to be able to filter the PARs)
  

foreach my $chromosome_slice (@{ $slice_a->fetch_all('toplevel', undef, 0, 1) }) {
  
  #Skip chr is we have specified just one to run
  if ($chr &&
      (uc($chr) ne uc($chromosome_slice->seq_region_name) )) {
    next;
  }

  #Get all black list regions for this slice
  my @excluded_regions = map $_->feature_Slice, 
    @{$chromosome_slice->get_all_MiscFeatures('encode_excluded')}; 
  

  #Add PARs for Y
  if ($chromosome_slice->seq_region_name eq 'Y') {
    #print "Adding the Y PAR regions to the blacklist\n";
    #We need to add them at the beginning so that these large regions will be tested first... 
    unshift @excluded_regions, $slice_a->fetch_by_region('chromosome','Y',10001,2649520);
    unshift @excluded_regions, $slice_a->fetch_by_region('chromosome','Y',59034050);
  }

  
  foreach my $region (@excluded_regions) {
    my @ids_to_remove;
    my %counts;
    my $par_cnt = 0;
    #could probably do with some global counts too.

    print "Fetching features from black list region:\t\t\t\t".$region->name."\n";
    my @set_features = @{$afa->fetch_all_by_Slice_FeatureSets($region, \@fsets)};
	
 
    #We need to maintain the  feature loop
    #if we want overlap stats and %blacklist

    #This assumes that anything projected is on a PAR
    #But could actually be projected from non top level slice!
    my $has_par = ($region->seq_region_name eq 'Y') ? 1 : 0;
	
    #Need to overhaul stats/logs here
 
    foreach my $feature (@set_features) {

      #Firstly avoid PAR regions which have been projected via the slice fetch code
      #$has_par is currently a whole chr test, rather than a test on the local region slice
      
      if ($has_par) {
        #Get the funcgen seq_region_id from the core Slice
        my $query_sr_id = $afa->get_seq_region_id_by_Slice($region);

        my $sql = 'select seq_region_id from annotated_feature where annotated_feature_id='.
          $feature->dbID;
        my ($db_sr_id) = $efgdba->dbc->db_handle->selectrow_array($sql);

        if ($query_sr_id != $db_sr_id){
          $par_cnt++;
          next;# if $query_sr_id != $db_sr_id;
        }
      }


      $global_counts{$feature->feature_set->name}{count} ++;
      $counts{$feature->feature_set->name}{count} ++;
      $counts{$feature->feature_set->name}{total_feature_length} += $feature->feature_Slice->length;


	  
      #This assumes that encode_excluded regions are non-overlapping with each other
      #There is an overlap
      # This has the disadvantage that overlap counts may be more than the features... 
      #$set_overlap_count++;  
      push @ids_to_remove, $feature->dbID();

      my $start = $region->start;
      my $end = $region->end;
	  
      if ($start < $feature->start) { 
        $start = $feature->start;  
      }
		
      if ($end > $feature->end) { 
        $end = $feature->end;  
      }
		
      #is this even useful?
      $counts{$feature->feature_set->name}{cumulative_region_coverage} += ($end-$start+1);

      #might be if we also had non_cumulative coverage

      #The same feature set can overlap several times with the same blacklist region...
      #Can probably update this or integrate with other log?
      $blacklist{$region->name}{$feature->feature_set->name} = 1;
    }
  

    foreach my $fset_name (keys %counts) {
      print FO "${fset_name}\t".join("\t", (map $counts{$fset_name}{$_}, @count_keys))."\n";
    }
    
    if(scalar(@ids_to_remove) == 0){
      print "No AnnotatedFeatures for blacklist region:\t".$region->name."\n";
    }
    elsif($remove){    
      &remove_features($region, \@ids_to_remove);
    }
    else{
      print "Skipping deletion of ".scalar(@ids_to_remove).
        " AnnotatedFeatures for region:\t".$region->name."\n";
    }
    
    if($par_cnt){
      print "Skipped $par_cnt features which were projected from a another seq_region (e.g. PAR)\n";
    }
  }
}


print FO "\n\nTotal counts\n";

foreach my $fset_name(keys %global_counts){

  print FO "${fset_name}\t".$global_counts{$fset_name}{count}."\n";
}


close FO;


open(FOB,">".$output.".blacklist");

foreach my $region (sort keys %blacklist){
  print FOB $region;
  
  foreach my $set (keys %{$blacklist{$region}}){ 
    print FOB "\t".$set; 
  }

  print FOB "\n";
}
close FOB;


sub remove_features{
  my ($region, $db_ids, $counts) = @_;

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

  return;
}




1;


