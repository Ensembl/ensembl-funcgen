#!/usr/bin/env perl

use strict;
use warnings;

use Test::More qw(no_plan); #no_plan required for skip usage

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

# switch on the debug prints
our $verbose = 0;
my $skip = 0;
my $okay = 0;

warn "FeatureType.t tests are incomplete and needs updating to use MultiTestDB\n";

ok(1, 'Startup test');#?

#my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
#my $db    = $multi->get_DBAdaptor( 'funcgen' );


#my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
#  (
#   -user    => 'XXX',
#   -host    => 'XXX',
#   -species => 'homo_sapiens', #Does this prevent alias loading?
#   -dbname  => 'XXX'
#  );

#debug( 'Test database instantiated' ); #Less verbose, but only get test names in line and in debug mode
ok( $db, 'DBAdaptor creation');# More verbose, but we get failed test name summary at end


my $dnadb = $db->dnadb;


#Add initial tests here for RegulatoryFeatureAdaptor
#e.g.
#list_dbIDs? only if we redefine this in funcgen in which case we need to rename
#the test accordingly? RegulatoryFeature_BaseAdaptor.t?
#other fetch methods
#and test attrs of known test data returned by these fetch methods
#see gene.t for examples

my $ftype_a = $db->get_FeatureTypeAdaptor;

ok($ftype_a, 'Got FeatureTypeAdaptor');

#TODO: {

#  todo_skip 'message here', 1;
#}  my $rfs     = $regf_a->fetch_all;


#Should get complete hashref
my $evidence_info = $ftype_a->get_regulatory_evidence_info;

ok( ($evidence_info && 
       (ref($evidence_info) eq 'HASH') ), 'Got regulatory evidence hash');

my @evidence_types = keys(%{$evidence_info});

#Simple content check
ok(scalar(@evidence_types) > 0, 'Regulatory evidence hash has keys');

#skip here if not


#my %regulatory_evidence_info = 
#  (
#   core => {
#            name      => 'Open chromatin & TFBS',
#            long_name => 'Open chromatin & Transcription factor binding sites',
#            label     => 'DNase1 & TFBS',
#            classes   => ['Transcription Factor', 'Transcription Factor Complex',   'Open Chromatin'],
#           },
#   
#   non_core => {
#                name      => 'Histones & polymerases',
#                long_name => 'Histone modifications & RNA polymerases',
#                label     => 'Hists & Pols',
#                classes   => ['Polymerase',  'Histone'],
#               }
#  );


#test $ftype_a->get_regulatory_evidence_info('BLART') warn and return nothing?
#test get_regulatory_evidence_classes or is this being removed?

SKIP: {

  if(! @evidence_types){
    skip 'Found no regulatory evidence config in FeatureTypeAdaptor', 5;
  }

  #Can assume the use of core and non_core
  #Can't avoid data test here as the web code uses this hash directly
  #Wrapper methods would avoid this
  #As it stands these data tests will have to be put into a HC
  #i.e. we should test for presence of name, long_name, label and classes
  #just implicitly check they match 


  ok(exists ${$evidence_info}{core},     'Regulatory evidence hash contains core key');
  ok(exists ${$evidence_info}{non_core}, 'Regulatory evidence hash contains non_core key');
  
  foreach my $evidence_type(@evidence_types){
    
    #Test we can use methods based on evidence_type and that they match the config
    #This requires bringing back an Ftype for each class and seeing whether it is core or non_core.
    my $evidence_type_info = $ftype_a->get_regulatory_evidence_info($evidence_type);

    #warn "evidence_type_info  $evidence_type $evidence_type_info\n".Data::Dumper::Dumper( $evidence_type_info)."\n";

    is_deeply($evidence_type_info, $evidence_info->{$evidence_type},
              "get_regulatory_evidence($evidence_type) matches full hash");


    #We have to check data here as we have no wrapper methods!!

    


    my $classes;

    if(exists $evidence_type_info->{classes}){
      $classes = $evidence_type_info->{classes};
    }
    
    my $got_classes = 0;
    
    if(defined $classes && 
       (ref($classes) eq 'ARRAY') && 
       (scalar(@$classes) > 0 ) ){
      $got_classes = 1;
    }

    ok($got_classes, "Got $evidence_type FeatureType classes");
    
  SKIP: {
      if(! $got_classes){
        skip "No $evidence_type FeatureType classes, cannot perform is_core test", 1;
      }
      
      foreach my $class(@$classes){
        my ($tmp_ftype) = @{$ftype_a->fetch_all_by_class($class)};
        $okay = 1;
        
        if($tmp_ftype->is_core_evidence &&
           ($evidence_type ne 'core') ){
          $okay = 0;
        }
        elsif( (! $tmp_ftype->is_core_evidence) && 
               ($evidence_type ne 'non_core') ){
          $okay = 0;
        }
        
        ok($okay, "$class FeatureType has correct is_core_evidence boolean");
      }
    }
  }
}

#Now we have check is_core_evidence works we can use it to validate return
#of fetch_all_by_evidence_type

my @ftypes;

#Test failure first
eval { @ftypes = $ftype_a->fetch_all_by_evidence_type('BLART') };
ok($@, 'fetch_all_by_evidence_type("BLART") failed');

@ftypes = @{$ftype_a->fetch_all_by_evidence_type('core')};

#Would really need to fetch all, and count those which are core 
#to validate we aren't missing anything?

$okay = 1;

foreach my $ftype(@ftypes){
  
  if($ftype->regulatory_evidence_type eq 'core'){
    $okay = 0;
    last;
  }
}

ok($okay, 'fetch_all_by_evidence_type("core") only returns core FeatureTypes');


@ftypes = @{$ftype_a->fetch_all_by_evidence_type('non_core')};

#Would really need to fetch all, and count those which are core 
#to validate we aren't missing anything?

$okay = 1;

foreach my $ftype(@ftypes){
  
  if($ftype->regulatory_evidence_type eq 'non_core'){
    $okay = 0;
  }

 SKIP: {
  
    if(! $okay ||
       $tested_evidence_methods ){
      skip '', 4;
    }
    
    
    #test regulatory_evidence_name/label/long_name here
    #and for core too?
    


  }


  last if ! $okay
  
}

ok($okay, 'fetch_all_by_evidence_type("non_core") only returns non_core FeatureTypes');








#done_testing();#was double printing, obviously called in Test::More DESTROY or something?
