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

if((scalar(@feature_sets)>0) && $all){ warn " Only specified feature sets will be considered. -all is being ignored."; }


#my @fsets = @{$fsa->fetch_all_by_type('annotated')};
my @fsets;
foreach my $feature_set (@feature_sets){ 
  my $fset = $fsa->fetch_by_name($feature_set);
  if($fset) { 
    push(@fsets, $fsa->fetch_by_name($feature_set)); 
  } else { warn "Could not find Feature Set $feature_set"; }
}
if(scalar(@fsets)==0){
  if($all){ 
    warn "All feature sets being used... this may take a while!";
    @fsets = @{$fsa->fetch_all()};
  } else { warn "No Feature Sets found. Use -feature_sets or -all"; exit 1; }
}
#my @fsets = @{$fsa->fetch_all_by_type('external')};
#my @fsets = @{$fsa->fetch_all()};
#print scalar(@fsets)."\n"; 

my %ids_to_remove;
my %blacklist;
open(FO,">".$output);
print FO "set_name\tset_count\tset_overlap_count\tset_length\tset_overlap_length\n";
foreach my $fset (@fsets){
  print $fset->name."\n";
  #print Dumper $fset;
  my $set_count = 0;
  my $set_overlap_count = 0;
  my $set_length = 0;
  my $set_overlap_length = 0;

  warn "Adding the Y PAR regions to the blacklist";
  my @all_regions = @{ $coredba->get_SliceAdaptor()->fetch_all('toplevel') };
  push @all_regions, $coredba->get_SliceAdaptor()->fetch_by_region('chromosome','Y',10001,2649521);
  push @all_regions, $coredba->get_SliceAdaptor()->fetch_by_region('chromosome','Y',59034050);


  foreach my $chromosome_slice (@{ $coredba->get_SliceAdaptor()->fetch_all('chromosome') }){
    #print "Chromosome\t".$chromosome_slice->seq_region_name()."\n";
    my @set_features = @{$fset->get_Features_by_Slice($chromosome_slice)};
    my @excluded_regions = @{$chromosome_slice->get_all_MiscFeatures('encode_excluded')}; 

    if($chromosome_slice->seq_region_name eq 'Y'){
      push @excluded_regions, $coredba->get_SliceAdaptor()->fetch_by_region('chromosome','Y',10001,2649521);
      push @excluded_regions, $coredba->get_SliceAdaptor()->fetch_by_region('chromosome','Y',59034050);
    }

    foreach my $feature (@set_features){
      $set_count++;
      $set_length += $feature->feature_Slice()->length();
      
      #This assumes that encode_excluded regions are non-overlapping with each other
      foreach my $region (@excluded_regions ) {
	if(($feature->start() >= $region->end()) || 
	   ($feature->end() <= $region->start()) ){ next; } else { 
	     #There is an overlap
	     # This has the disadvantage that overlap counts may be more than the features... 
	     $set_overlap_count++;

	     #potential memory hog... there could be many features... but should be ok.
	     $ids_to_remove{$feature->dbID()} = 1;

	     my $start = $region->start();
	     my $end = $region->end();
	     if($start<$feature->start()){ $start = $feature->start();  }
	     if($end>$feature->end()){ $end = $feature->end();  }
	     $set_overlap_length += ($end-$start+1);

	     #The same feature set can overlap several times with the same blacklist region...
	     $blacklist{$chromosome_slice->seq_region_name()."_".$region->start()."_".$region->end()}{$fset->name} = 1;

	   }	
      }
    }
  }
  print FO $fset->name."\t".$set_count."\t".$set_overlap_count."\t".$set_length."\t".$set_overlap_length."\n";
}
close FO;


open(FOB,">".$output.".blacklist");
foreach my $region (sort keys %blacklist){
  print FOB $region;
  foreach my $set (keys %{$blacklist{$region}}){ print FOB "\t".$set; }
  print FOB "\n";
}
close FOB;

if($remove){
  #TODO Batch remove features!
  foreach my $fID (keys %ids_to_remove){
    eval{ 
      my $sql = "DELETE FROM annotated_feature WHERE annotated_feature_id=".$fID.";";
      $efgdba->dbc->do($sql);		 
    };
    if($@) { warn $@->getErrorMessage(); }
  }
} 

exit 0;


