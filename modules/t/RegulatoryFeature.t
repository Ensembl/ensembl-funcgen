# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2017] EMBL-European Bioinformatics Institute
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
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

# switch on the debug prints
our $verbose = 0;
my $skip = 0;

ok(1, 'Startup test');#?

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db    = $multi->get_DBAdaptor( 'funcgen' );



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

my $regf_a = $db->get_RegulatoryFeatureAdaptor;

TODO: {

  todo_skip 'RegulatoryFeatureAdaptor->fetch_all needs small test DB data set', 1;
  my $rfs     = $regf_a->fetch_all;
  my $success = 1;

  if(ref($rfs) ne 'ARRAY'){
    $success = 0;
  }
  else{
   foreach my $rf(@$rfs){
     if(! (defined($rf) && 
           ref($rf) eq 'Bio::EnsEMBL::Funcgen::RegulatoryFeature') ){
       $success = 0;
     }
   }
  }
  
  okay($success, 'RegulatoryFeatureAdaptor::fetch_all');
}



### BOUND/PAR TESTS ###


# TO DO
# 1 Test with -ve query slice
# 2 Set up test DB with known data
# 3 Test a feature where the attrs do not exceed the core seq_region loci
# 4 Add get_underlying_structure tests as this has changed

# now create/fetch a new RegulatoryFeature on X PAR
my $slice_a = $dnadb->get_SliceAdaptor;
#my $x_slice = $slice_a->fetch_by_region('chromosome', 'X', 68800, 72000, -1);
my $x_start = 332794;
my $x_end   = 380000;
my $x_slice = $slice_a->fetch_by_region('chromosome', 'X', $x_start, $x_end);


#Use H1ESC as this seems to have the most bounds in this region
#We still have issues around overlapping data here!
#How were these not caught/fitlered?
#my ($rf) = @{$fset->get_Features_by_Slice($x_slice)};
my $rf;
my $fset = $db->get_FeatureSetAdaptor->fetch_by_name('RegulatoryFeatures:H1ESC');

foreach my $regf(@{$fset->get_Features_by_Slice($x_slice)}){
  #warn $regf->stable_id.' '.$regf->bound_start.' - '.$regf->start.':'.$regf->end.' - '.$regf->bound_end."\n";

  if( ($regf->bound_seq_region_start != $regf->seq_region_start) &&
      ($regf->bound_seq_region_end   != $regf->seq_region_end) ){
    $rf = $regf;
    last;
  }
}

#Sanity check we actually have bounds
#Could have this as this as the SKIP condition 
#But this is also needs to be a test
$skip = 1 if ! defined $rf;
ok( $rf, 'Found RegulatoryFeature with bounds');

SKIP: {
  #Could have debug here, but this only print if we skip anyway
  if( $skip ){
    #skip 'RegulatoryFeature('.$rf->dbID.") bounds do not differ from seq_region loci:\t".
    #$rf->bound_seq_region_start.' - '. $rf->seq_region_start.' -- '.
    #    $rf->seq_region_end.' - '. $rf->bound_seq_region_end;

    skip('Could not identify RegulatoryFeature with start and end bound', 13);
  }

   
  # Need to revert funcgen BaseFeatureAdaptor work around first?

  #Let's test they are < or > start/end

  ok($rf->bound_seq_region_start < $rf->seq_region_start, 
     'bound sr start < sr start');
  
  ok($rf->bound_start < $rf->start,
       'bound local start > start');

  ok($rf->bound_seq_region_end > $rf->seq_region_end, 
     'bound sr end > sr end');

  ok($rf->bound_end > $rf->end,
     'bound local end > end');


 #TODO:{
  #local $TODO = 'Implement bound_start/end_length methods';
    #Doesn't work unless the code exists!
    #todo_skip 'bound_start/end_length methods not yet implemented', 4;
  #}


  my @struc = @{$rf->get_underlying_structure};

  ok( ($rf->bound_start == $struc[0]) &&
      ($rf->start       == $struc[1]) &&
      ($rf->end         == $struc[-2]) &&
      ($rf->bound_end   == $struc[-1]),
      'underying_structure bound/start/ends match');

  &test_bound_length_start_end($rf); 
  

  # Now get respective Y PAR feature
 
  

  my $y_slice = $slice_a->fetch_by_region('chromosome', 'Y', ($x_start - 50000), ($x_end - 30000) );
  #use x seq_region_end as this is always ~40kb > Y seq_region_end
  my $y_rf;

  foreach my $reg_feat(@{$fset->get_Features_by_Slice($y_slice)}){
    
    if($reg_feat->dbID == $rf->dbID){
      $y_rf = $reg_feat;
      last;
    }
  }

  # Alternative would be to project to Y PAR
  # via Slice::project_to_slice or Feature::transfer/project_to_slice
  # But these do not yet support HAP/PAR mapping (only assembly/coord_system)

  ok($y_rf, 'Corresponding Y PAR RegulatoryFeature fetched');


  SKIP: {

    if(! defined $y_rf){
      skip('Skipping Y par tests as failed to fetch corresponding Y PAR feature', 3);
    }

    ok(($rf->seq_region_start != $y_rf->seq_region_start) &&
       ($rf->seq_region_end != $y_rf->seq_region_end),
       'Y RegulatoryFeature projection seq_region loci do not match X');
  
    #These are proxy tests until we know the exact data/values
    ok($y_rf->bound_start != $rf->bound_start, 'projected bound_start changed');
    ok($y_rf->bound_end   != $rf->bound_end, 'projected bound_start changed');
  }
}

$skip = 0; 


sub test_bound_length_start_end{
  my $rf = shift;
  
  #Using a -ve slice causes these to fail due to 
  #local boundl lengths

  #This is essentially recreating the API calc!
  #Which is not what we want to do.
  #Should test known data

  ok($rf->bound_start_length ==
     ($rf->seq_region_start - $rf->bound_seq_region_start),
     'bound_start_length seq_region calc');
  
  ok($rf->bound_start_length ==
     ($rf->start - $rf->bound_start),
     'bound_start_length local calc');
  
  ok($rf->bound_end_length ==
     ($rf->bound_seq_region_end - $rf->seq_region_end),
     'bound_end_length seq_region calc');
  
  ok($rf->bound_end_length ==
     ($rf->bound_end - $rf->end),
     'bound_end_length local calc');

  return;
}


# Old vs New Build tests
# These are for methods which return different data
# from different build versions e.g. 
# summary_as_hash, has_evidence, cell_type_count & is_projected
# So these need to be done with two builds
# Here we will use 2 dbs for convinience
# but this needs changing to use different feature sets
# once we start using the test DB.


#Assuming that the last $rf we saw was from a new build version


test_summary_as_hash($rf, 'New');

my $mdb = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
  (-user    => 'XXX',
   -host    => 'XXX',
   -species => 'mus_musculus', #Does this prevent alias loading?
   -dbname  => 'mus_musculus_funcgen_78_38'
  );

my $sid = 'ENSMUSR00000233228';
my $mrf = $mdb->get_RegulatoryFeatureAdaptor->fetch_by_stable_id($sid);

SKIP: {

  if(! defined $mrf){
    skip("Skipping Old build tests as failed to retrieve $sid", 1);
  }

  test_summary_as_hash($mrf, 'Old');
}


sub test_summary_as_hash{
  my $rf         = shift;
  my $build_type = shift;

  my $hash_summary = $rf->summary_as_hash;

  #my @hash_array = %$hash_summary;
  #Checking for any undefs translated to missing elements
  #ok(! (scalar(@hash_array) % 2), 'summary_as_hash returns even sized list');
  #This will always return an even sized list as perl will have simply shifted things
  #up into the apparent void and appended and undef


  #Now let's check the ones we know about and for any unknown ones?
  #although it may be better to omit undef kv pairs
  #in terms of REST performance, we want to reduce the amount of data return
  #and let the calling API/code handle/translate the ommissions as undefs.

  my %summary_keys = 
   (ID                => undef,
    cell_type         => undef,
    bound_start       => undef, 
    bound_end         => undef, 
    start             => undef, 
    end               => undef, 
    strand            => undef, 
    seq_region_name   => undef, 
    activity_evidence => undef, 
    description       => undef, 
    feature_type      => undef,
    projected         => undef,
    cell_type_count   => undef);

  my $hash_valid = 1; 

  foreach my $key(keys %$hash_summary){

    if( ! exists $summary_keys{$key}){
      $hash_valid = 0;
      last;
    }
  }

  my $invalid_keys = ($hash_valid) ? '' : 
   join(' ', keys %$hash_summary);

  ok($hash_valid, "Summary hash invalid keys:\t$invalid_keys");

}

#done_testing();#was double printing, obviously called in Test::More DESTROY or something?
