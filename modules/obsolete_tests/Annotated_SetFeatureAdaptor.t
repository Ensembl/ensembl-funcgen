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

warn "Annotated_SetFeatureAdaptor.t tests are incomplete and needs updating to use MultiTestDB\n";


ok(1, 'Startup test');#?

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db    = $multi->get_DBAdaptor( 'funcgen' );



#debug( 'Test database instantiated' ); #Less verbose, but only get test names in line and in debug mode
ok( $db, 'DBAdaptor creation');# More verbose, but we get failed test name summary at end

my $af_a    = $db->get_AnnotatedFeatureAdaptor;
ok($af_a, 'Got AnnotatedFeatureAdaptor');

my $fset_a  = $db->get_FeatureSetAdaptor;
my $slice_a = $db->dnadb->get_SliceAdaptor;
my $slice   = $slice_a->fetch_by_region('chromosome', 1, 5000000, 10000000);
my $ftype_name = 'CTCF';
my $ftype   = $db->get_FeatureTypeAdaptor->fetch_by_name($ftype_name);

#TODO update skip counts

SKIP: {

  if(! $ftype){
    skip "Found no $ftype_name FeatureType, please amend this test", 8;
  }

  #SetFeatureAdaptor::fetch_all_by_Slice_FeatureType (optional logic_name)

  my $feats = $af_a->fetch_all_by_Slice_FeatureType($slice, $ftype);
  ok(ref($feats) eq 'ARRAY', "fetch_all_by_Slice_FeatureType - Fetched $ftype_name features");
  
  SKIP: {
    
    if(! scalar(@$feats)){
      skip "Found no $ftype_name features found, please amend this test", 7;      
    }
    
    my $all_feats = scalar(@$feats);
    #warn "all feates $all_feats";
    my $feat_cnt    = 0;
    
    #TODO add test for feature count here
    
    #Could do this without knowledge of DB contents by getting all CTCF fsets
    #and their logic names, then querying for each, and making sure they add up.
    
    my @fsets = @{$fset_a->fetch_all_by_FeatureType($ftype)};
    my %lnames;
    map {$lnames{$_->analysis->logic_name} = undef} @fsets; 
    foreach my $lname(keys %lnames){
      $feat_cnt += scalar(@{$af_a->fetch_all_by_Slice_FeatureType($slice, $ftype, $lname)});
      #warn "feat_cnt now $feat_cnt after $lname query";
    }
    
    
    ok($all_feats == $feat_cnt, 'fetch_all_by_Slice_FeatureType - Got all features using logic_names');
    
    $feats = $af_a->fetch_all_by_Slice_FeatureType($slice, $ftype, 'NONE EXISTANT LOGIC NAME');   
    ok(scalar(@$feats) == 0, 
       'fetch_all_by_Slice_FeatureType - None existant logic_name returned no features');
    
    #SetFeatureAdaptor::fetch_all_by_FeatureSets (optional logic_name);
    my $fset   = $fsets[0];
    my $lname  = $fset->analysis->logic_name;   
    $feats     = $af_a->fetch_all_by_FeatureSets([$fset]);
    $all_feats = scalar(@$feats);
    ok($all_feats > 0, 'fetch_all_by_FeatureSets - Fetched FeatureSet features');
    
    $feats     = $af_a->fetch_all_by_FeatureSets([$fset], $lname);
    ok($all_feats == scalar(@$feats), 
       'fetch_all_by_FeatureSets - Fetched FeatureSet features with logic_name');
    
    #No need for NON EXISTANT LOGIC NAME here as we have already tested it
    #In fact, isn't that a BaseAdaptor test?
    #We just want to test that logic_names are passed on correctly
    
    #SetFeatureAdaptor::fetch_all_by_Slice_FeatureSets (optional logic_name);
    $feats     = $af_a->fetch_all_by_Slice_FeatureSets($slice, [$fset]);
    $all_feats = scalar(@$feats);
    ok($all_feats > 0, 'fetch_all_by_Slice_FeatureSets - Fetched FeatureSet Slice features');
    
    $feats     = $af_a->fetch_all_by_Slice_FeatureSets($slice, [$fset], $lname);
    ok($all_feats == scalar(@$feats), 
       'fetch_all_by_Slice_FeatureSets - Fetched FeatureSet Slice features with logic_name');

    #SetFeatureAdaptor::fetch_Iterator_by_Slice_FeatureSets (optional logic_name and chunk length);
    my $iter = $af_a->fetch_Iterator_by_Slice_FeatureSets($slice, [$fset]);
    ok(check_ref($iter, 'Bio::EnsEMBL::Utils::Iterator'),
       'fetch_Iterator_by_Slice_FeatureSets - Got Iterator');
    
    
    #what else do we need to check here wrt iterator method

  }
}
    


#done_testing();#was double printing, obviously called in Test::More DESTROY or something?
