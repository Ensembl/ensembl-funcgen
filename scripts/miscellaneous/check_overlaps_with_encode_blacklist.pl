#!/software/bin/perl

=head1 NAME

check_overlaps_with_blacklist.pl - For all the human annotated features, checks how they overlap with a set of locations defined as being problematic by the ENCODE project (so called blacklist)

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

=cut

use warnings;
use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Data::Dumper;
use Getopt::Long;
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
my $remove;
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
            "remove"        => \$remove,
            "feature_sets=s{,}"  => \@feature_sets);

if(scalar(@feature_sets)==0){ throw("Must specify feature set(s) name with -feature_sets option"); exit 1; }

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

my $fsa = $efgdba->get_FeatureSetAdaptor();
my $afa = $efgdba->get_AnnotatedFeatureAdaptor();
#my $cta = $efgdba->get_CellTypeAdaptor();
#my $K562 = $cta->fetch_by_name("K562");
#This doesn't seem to be working very well... TODO: check why...
#my @fsets = $fsa->fetch_all_by_CellType($K562);

#my @fsets = @{$fsa->fetch_all_by_type('annotated')};
my @fsets;
foreach my $feature_set (@feature_sets){ 
  my $fset = $fsa->fetch_by_name($feature_set);
  if($fset) { 
    push(@fsets, $fsa->fetch_by_name($feature_set)); 
  } else { warn "Could not find Feature Set $feature_set"; }
}
if(scalar(@fsets)==0){ warn "No Feature Sets found"; exit 1; }
#my @fsets = @{$fsa->fetch_all_by_type('external')};
#my @fsets = @{$fsa->fetch_all()};
#print scalar(@fsets)."\n"; 

open(FO,">".$output);
foreach my $fset (@fsets){
  print $fset->name."\n";
  #print Dumper $fset;
  my $set_count = 0;
  my $set_overlap_count = 0;
  my $set_length = 0;
  my $set_overlap_length = 0;
  foreach my $chromosome_slice (@{ $coredba->get_SliceAdaptor()->fetch_all('chromosome') }){
    #print "Chromosome\t".$chromosome_slice->seq_region_name()."\n";
    my @set_features = @{$fset->get_Features_by_Slice($chromosome_slice)};
    my @excluded_regions = @{$chromosome_slice->get_all_MiscFeatures('encode_excluded')}; 

    foreach my $feature (@set_features){
      $set_count++;
      $set_length += $feature->feature_Slice()->length();
      
      #This assumes the encode_excluded regions are non-overlapping with each other
      foreach my $region (@excluded_regions ) {
	if(($feature->start() >= $region->end()) || 
	   ($feature->end() <= $region->start()) ){ next; } else { 
	     #There is an overlap
	     # This has the disadvantage that overlap counts may be more than the features... 
	     $set_overlap_count++;

	     my $start = $region->start();
	     my $end = $region->end();
	     if($start<$feature->start()){ $start = $feature->start();  }
	     if($end>$feature->end()){ $end = $feature->end();  }
	     $set_overlap_length += ($end-$start+1);

             if($remove){
	       eval{ 
		 #print "Removing Feature with DB ID ".$feature->dbID()."\n";
		 my $sql = "DELETE FROM annotated_feature WHERE annotated_feature_id=".$feature->dbID().";";
		 $efgdba->dbc->do($sql);		 
	       };
	       if($@) { warn $@->getErrorMessage(); }
	       last; #pass to the next feature... 
	       #overlap counts will be more accurate but overlap length may not be accurate in this case...
	     } 

	   }	
      }
    }
  }
  print FO $fset->name."\t".$set_count."\t".$set_overlap_count."\t".$set_length."\t".$set_overlap_length."\n";
}
close FO;

exit 0;


