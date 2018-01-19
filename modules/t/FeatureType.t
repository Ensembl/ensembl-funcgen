# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.



use strict;
use warnings;

use Test::More qw(no_plan); #no_plan required for skip usage

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Utils::Scalar qw( check_ref );  
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;


# switch on the debug prints
our $verbose = 0;
my $skip = 0;
my $okay = 0;

ok(1, 'Startup test');#?

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db    = $multi->get_DBAdaptor( 'funcgen' );
my $dnadb = $db->dnadb;


#Add initial tests here for RegulatoryFeatureAdaptor
#e.g.
#list_dbIDs? only if we redefine this in funcgen in which case we need to rename
#the test accordingly? RegulatoryFeature_BaseAdaptor.t?
#other fetch methods
#and test attrs of known test data returned by these fetch methods
#see gene.t for examples

my $ftype_a = $db->get_FeatureTypeAdaptor;

isa_ok($ftype_a, 'Bio::EnsEMBL::Funcgen::DBSQL::FeatureTypeAdaptor', 'Got FeatureTypeAdaptor');

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
  
  if($ftype->regulatory_evidence_type ne 'core'){
    $okay = 0;
    last;
  }
}

ok($okay, 'fetch_all_by_evidence_type("core") only returns core FeatureTypes');


@ftypes = @{$ftype_a->fetch_all_by_evidence_type('non_core')};

#Would really need to fetch all, and count those which are core 
#to validate we aren't missing anything?

$okay = 1;
my $tested_evidence_methods = 0;

foreach my $ftype(@ftypes){
  
  if($ftype->regulatory_evidence_type ne 'non_core'){
    $okay = 0;
  }

  #Do we need to do this for core too?
  if(! $tested_evidence_methods){
    $tested_evidence_methods = 1;
    

  SKIP: {
      #This should only skip if the evidence_type is not defined?

      if(! defined $ftype->regulatory_evidence_type){
        skip 'Cannot test regulatory evidence methods with undefined evidence_type', 3;
      }
    
      #add string test here?
      my $evidence = $ftype->evidence_type_name;
      ok(defined $evidence &&
         (ref(\$evidence) eq 'SCALAR'), 'evidence_type_name defined');

      $evidence = $ftype->evidence_type_long_name;
      ok(defined $evidence &&
         (ref(\$evidence) eq 'SCALAR'), 'evidence_type_long_name defined');
      
      $evidence = $ftype->evidence_type_label;
      ok(defined $evidence &&
         (ref(\$evidence) eq 'SCALAR'), 'evidence_type_label defined');
    }
  }

  last if ! $okay;
}

ok($okay, 'fetch_all_by_evidence_type("non_core") only returns non_core FeatureTypes');


#Test redundant name fetching. 
my $ftype;

eval{ $ftype = $ftype_a->fetch_by_name('Predicted Transcribed Region') };
ok($@, 'fetch_by_name caught multiple FeatureTypes returned to a scalar context');

eval{ @ftypes = $ftype_a->fetch_by_name('Predicted Transcribed Region') };
my $error = $@;

SKIP: {
      #This should only skip if the evidence_type is not defined?

      if(scalar(@ftypes) < 2){
        skip 'Cannot test fetch_by_name returns array as there are <1 Predicted Transcribed Region feature types stored', 1;
      }

      ok((! $error) && check_ref($ftypes[0], 'Bio::EnsEMBL::Funcgen::FeatureType'), 'fetch_by_name returns multiple FeatureTypes in a scalar context');
}

my $ftypes;
eval{ $ftypes = $ftype_a->fetch_all_by_name('Predicted Transcribed Region') };
$error = $@;

#We should load these here, instead of depending 
#on pre-stored Predicted Transcribed Region entries.
# ok((! $@) && 
#    check_ref($ftypes, 'ARRAY') &&
#    scalar(@$ftypes) == 24      &&
#    check_ref($ftypes->[0], 'Bio::EnsEMBL::Funcgen::FeatureType'), 'fetch_all_by_name returns expected number(expected 24, got '.scalar(@$ftypes).') of FeatureTypes');










