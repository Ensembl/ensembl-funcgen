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
use Bio::EnsEMBL::Utils::Scalar            qw( check_ref );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( are_valid );
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

# switch on the debug prints
our $verbose = 0;
my $skip = 0;
my $okay = 0;

warn "Annotated_SetFeatureAdaptor.t tests are incomplete and needs updating to use MultiTestDB\n";


ok(1, 'Startup test');#?

#my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
#my $db    = $multi->get_DBAdaptor( 'funcgen' );


my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
  (
   -user    => 'XXX',
   -host    => 'XXX',
   -species => 'homo_sapiens',
   -dbname  => 'homo_sapiens_funcgen_73_37'
  );


#debug( 'Test database instantiated' ); #Less verbose, but only get test names in line and in debug mode
ok( $db, 'DBAdaptor creation');# More verbose, but we get failed test name summary at end


SKIP : {
  
  if(! $db){
    skip 'Failed to create DBAdaptor, skipping all tests', 0;
  }

  my $iset_a = $db->get_InputSetAdaptor;
  ok($iset_a, 'Got InputSetAdaptor');
 

  #TODO update skip counts
  #Add test for optional status arg to fetch_all_by_Feature/CellType

  SKIP: {
    if(! defined $iset_a){
      skip 'Failed to get_InputSetAdaptor', 0;  
    }
  
    #Get all isets via fetch_all and cache, then we can test numbers returned from
    #fetch methods without knowing the contents of the DB apiori
    #This is still not perfect, but provides some level of rigour
    #until we have a known test DB
    my ($table, $table_syn) = @{$iset_a->_main_table};
    
    ok((($table eq 'input_set') &&
        ($table_syn eq 'inp')),
      'BaseAdaptor::_main_table gives correct name and synonym');
    
    my $dbIDs = $iset_a->_list_dbIDs;
    ok(check_ref($dbIDs, 'ARRAY'), 'BaseAdaptor::_list_dbIDs gives Arrayref');
       
    my $isets = $iset_a->fetch_all;
    eval { are_valid('Bio::EnsEMBL::Funcgen::InputSet', $isets) };
    my $failed = ($@) ? 1 : 0;
    ok( ((! $failed) &&
          scalar(@$isets) == scalar(@$dbIDs)),
       'BaseAdaptor::fetch_all InputSets matches number from _list_dbIDs');
    
    SKIP: {    
      if($failed){
        skip 'Found no InputSets found, please amend/fix this test', 5;
      }

      my %iset_cache;
      
      foreach my $iset(@$isets){
        $iset_cache{feature_type}{$iset->feature_type->name} ||= [];
        push @{$iset_cache{feature_type}{$iset->feature_type->name}}, $iset;
        
        $iset_cache{cell_type}{$iset->cell_type->name} ||= [];
        push @{$iset_cache{cell_type}{$iset->cell_type->name}}, $iset;
        
        #This is currently 1:1 but may change in future
        $iset_cache{experiment}{$iset->get_Experiment->name} ||= [];
        push @{$iset_cache{experiment}{$iset->get_Experiment->name}}, $iset;
      }


      #Absent argument is actually a BaseAdaptor test, as this is actually a test on
      #compose_constraint_query and the individual constrain methods
      #These can be done directly on the methods
      #not need for an actual query
      #so long as the fetch methods here return the expected results

      #InputSetAdaptor::fetch_all_by_FeatureType
      my ($ftype_name) = ( keys(%{$iset_cache{feature_type}}) ); #Take first one
      my $ftype        = $iset_cache{feature_type}{$ftype_name}->[0]->feature_type;
      $isets           = $iset_a->fetch_all_by_FeatureType($ftype);
      eval { are_valid('Bio::EnsEMBL::Funcgen::InputSet', $isets) };
      ok(((! $@) && 
          (scalar(@$isets) eq scalar(@{$iset_cache{feature_type}{$ftype_name}})) ), 
         "fetch_all_by_FeatureType - Fetched correct amount of InputSets");
      
      #(Input)SetAdaptor::_constrain_feature_types
      eval{ $iset_a->_constrain_feature_types([]) };
      $failed = ($@) ? 1 : 0;
      ok($failed, '_constrain_feature_types - Empty arg Arrayref');
      
      eval{ $iset_a->_constrain_feature_types() };
      $failed = ($@) ? 1 : 0;
      ok($failed, '_constrain_feature_types - No args');  
       
      
      #InputSetAdaptor::fetch_all_by_CellType
      my ($ctype_name) = ( keys(%{$iset_cache{cell_type}}) ); #Take first one
      my $ctype        = $iset_cache{cell_type}{$ctype_name}->[0]->cell_type;
      $isets           = $iset_a->fetch_all_by_CellType($ctype);
      eval { are_valid('Bio::EnsEMBL::Funcgen::InputSet', $isets) };
      ok(((! $@) && 
          (scalar(@$isets) eq scalar(@{$iset_cache{cell_type}{$ctype_name}})) ), 
         "fetch_all_by_CellType - Fetched correct amount of InputSets");
      
      #(Input)SetAdaptor::_constrain_cell_types
      eval{ $iset_a->_constrain_cell_types([]) };
      $failed = ($@) ? 1 : 0;
      ok($failed, '_constrain_cell_types - Empty arg Arrayref');
      
      eval{ $iset_a->_constrain_cell_types() };
      $failed = ($@) ? 1 : 0;
      ok($failed, '_constrain_cell_types - No args');
        
      
      #InputSetAdaptor::fetch_all_by_Experiment
      my ($exp_name) = ( keys(%{$iset_cache{experiment}}) ); #Take first one
      my $exp        = $iset_cache{experiment}{$exp_name}->[0]->get_Experiment;
      $isets         = $iset_a->fetch_all_by_Experiment($exp);
      eval { are_valid('Bio::EnsEMBL::Funcgen::InputSet', $isets) };
      ok(((! $@) && 
          (scalar(@$isets) eq scalar(@{$iset_cache{experiment}{$exp->name}})) ), 
         "fetch_all_by_Experiment - Fetched correct amount of InputSets");
       
      #(Input)SetAdaptor::_constrain_experiment
      eval{ $iset_a->_constrain_experiments([]) };
      $failed = ($@) ? 1 : 0;
      ok($failed, '_constrain_experiments - Empty arg Arrayref');
      
      eval{ $iset_a->_constrain_experiments() };
      $failed = ($@) ? 1 : 0;
      ok($failed, '_constrain_experiments - No args');
      
          
      #InputSetAdaptor::fetch_by_name
      eval {$iset_a->fetch_by_name()};
      $failed = ($@) ? 1 : 0;
      ok($failed, 'fetch_by_name - Caught no name argument');
            
      my $iset1 = $isets->[0];
      my $iset2 = $iset_a->fetch_by_name($iset1->name);
      ok((check_ref($iset2, 'Bio::EnsEMBL::Funcgen::InputSet') && 
          ($iset1->dbID == $iset2->dbID)),
         "fetch_by_name - Fetched correct InputSet");
      
      
      #Not testing _constraint_format yet as this is likely going to be removed
           
      #Add in analysis stuff here after patch
      
      
      
      #BaseAdaptor tests
      #(Input)SetAdaptor::_constrain_experiment
      eval{ $iset_a->_constrain_states([]) };
      $failed = ($@) ? 1 : 0;
      ok($failed, 'BaseAdaptor::_constrain_states - Empty arg Arrayref');
      
      eval{ $iset_a->_constrain_states() };
      $failed = ($@) ? 1 : 0;
      ok($failed, 'BaseAdaptor::_constrain_states - No args');
      
      
    }
  }
}
    


#done_testing();#was double printing, obviously called in Test::More DESTROY or something?
