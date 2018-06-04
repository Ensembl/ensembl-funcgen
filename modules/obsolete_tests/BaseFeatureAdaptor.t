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

ok(1, 'Startup test');

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $db    = $multi->get_DBAdaptor( 'funcgen' );


#debug( 'Test database instantiated' ); #Less verbose, but only get test names in line and in debug mode
ok( $db, 'DBAdaptor creation');# More verbose, but we get failed test name summary at end


#Use ProbeFeatureAdaptor as MotifFeatureADaptor may become a SetFeatureAdaptor
#if we host below threshold motif mapping sets

my $pf_a    = $db->get_ProbeFeatureAdaptor;
ok($pf_a, 'Got ProbeFeatureAdaptor');
my $slice_a = $db->dnadb->get_SliceAdaptor;
my $slice   = $slice_a->fetch_by_region('chromosome', 1, 9000000, 10000000);

#TODO update skip counts
my $lname     = 'AFFY_UTR_ProbeAlign';
my $analysis  = $db->get_AnalysisAdaptor->fetch_by_logic_name($lname);


SKIP: {

  if(! $analysis){
    skip "Found no $lname analysis, please amend this test or test DB", 8;
  }

  my $feats = $pf_a->fetch_all_by_Slice_constraint($slice, 'pf.analysis_id='.$analysis->dbID);
  ok(check_ref($feats, 'ARRAY'), 
   'fetch_all_by_Slice_constraint - Constraint string returns Arrayref');
  

  my $param_feats = 
   $pf_a->fetch_all_by_Slice_constraint($slice, {constraints => {logic_names =>[$lname]}});
  ok(check_ref($feats, 'ARRAY'), 
   'fetch_all_by_Slice_constraint - Constraint param returns Arrayref');  


  SKIP: {
    
    if(! (check_ref($feats, 'ARRAY') && check_ref($param_feats, 'ARRAY'))){
      skip "Cannot compare string/param contrained feature counts due to previous failure", 1;      
    }
    
    my $feat_cnt  = scalar(@$feats);
    my $num_feats = scalar(@$param_feats);
    ok($feat_cnt == $num_feats,
     'fetch_all_by_Slice_constraint - '.
      "String/param constrained queries return same amount of features ($feat_cnt vs $num_feats)");
  }

  my $lname_feats = 
   $pf_a->fetch_all_by_Slice_constraint($slice, undef, $lname);
  ok(check_ref($lname_feats, 'ARRAY'), 
   'fetch_all_by_Slice_constraint - Constraint logic_name returns Arrayref');  

  SKIP: {
    
    if(! (check_ref($feats, 'ARRAY') && check_ref($lname_feats, 'ARRAY'))){
      skip "Cannot compare string/logic_name contrained feature counts due to previous failure", 1;      
    }

    my $feat_cnt  = scalar(@$feats);
    my $num_feats = scalar(@$lname_feats);
    ok($feat_cnt == $num_feats,
     'fetch_all_by_Slice_constraint - '.
      "String/logic_name constrained queries return same amount of features ($feat_cnt vs $num_feats)");
  }
}
    


#done_testing();#was double printing, obviously called in Test::More DESTROY or something?
