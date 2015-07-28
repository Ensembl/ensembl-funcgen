# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#!usr/bin/env perl

use strict;
use warnings;

use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Test::More;
use Test::Exception;               # throws_ok
use Bio::EnsEMBL::Test::TestUtils  qw( test_getter_setter debug );
use Data::Dumper                   qw( Dumper );
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Test::MultiTestDB;

ok(1, 'Start up');

# switch on the debug prints
our $verbose = 0;

my $multi  = Bio::EnsEMBL::Test::MultiTestDB->new();
my $efgdba = $multi->get_DBAdaptor("funcgen");

# This test uses the following comment conventions
# START method_name testing
# COMPLETED method_name testing



my $dsa = $efgdba->get_adaptor("dataset");

#Just grab a few data sets to work with
#DISPLAYABLE ensure we should have an input_set defined
#$dsa->fetch_all;
my $dset   = $dsa->fetch_by_name('RegulatoryFeatures:MultiCell');
my $dset_2 = $dsa->fetch_by_name('RegulatoryFeatures:NHEK');

if(! (defined $dset && defined $dset_2)){
  throw('Failed to fetch 2 DataSets to test, please update DataSet.t');  
}


# START testing compare_to

my $diffs = '';
my %diffs = %{$dset->compare_to($dset)};
$diffs = "\n".Dumper(\%diffs) if %diffs;
ok(! %diffs, 'DataSet::compare_to self default no diffs'.$diffs);

# redefine methods
%diffs = %{$dset->compare_to($dset, undef,
  [qw(name get_all_states)] ,
  [qw(product_FeatureSet get_supporting_sets)]  
  )};
$diffs = "\n".Dumper(\%diffs) if %diffs;
ok(! %diffs, 'DataSet::compare_to self redefined methods no diffs'.$diffs);

eval{#Invalid redefined scalar method
 $dset->compare_to($dset, undef, 
                   [qw(product_FeatureSet get_supporting_sets)] ); 
};
ok($@, 'DataSet::compare_to scalar methods redefined (invalid)');


eval{#Invalid redefined object method
 $dset->compare_to($dset, undef, undef,
                   [qw(name get_all_states)] ); 
};
ok($@, 'DataSet::compare_to object methods redefined (invalid)');

# COMPLETED testing compare_to

## START testing reset_relational_attributes                                             

#use eq directly on object refs instead of compare_to


#clone ResultSet, so we can change some attrs
my $clone_dset = bless({%{$dset}}, ref($dset));
my $alt_fset   = $dset_2->product_FeatureSet;
my $alt_ssets  = $dset_2->get_supporting_sets;

my %relational_params = 
  (
   -feature_set     => $alt_fset,
   -supporting_sets => $alt_ssets,
  );


$clone_dset->reset_relational_attributes(\%relational_params, 'no_db_reset');


#eq comparisons of obj in scalar context compares mem refs
ok( ($clone_dset->product_FeatureSet eq $alt_fset),
   'DataSet::reset_relational_attributes reset FeatureSet');

#can't compare arrayrefs directly, so use compare_to here
%diffs = %{$clone_dset->compare_to($dset_2, undef, [],
                                      ['get_supporting_sets'])};
$diffs = "\n".Dumper(\%diffs) if %diffs;
ok(! $diffs,
   'DataSet::reset_relational_attributes reset supporting sets'.$diffs);


ok((defined $clone_dset->dbID && 
   defined $clone_dset->adaptor), 
   'DataSet::reset_relational_attributes no_db_reset');
    
$clone_dset->reset_relational_attributes(\%relational_params);
ok(! (defined $clone_dset->dbID || 
      defined $clone_dset->adaptor), 
      'DataSet::reset_relational_attributes with dbID/adaptor reset');



eval { $clone_dset->reset_relational_attributes(
        {  
         -feature_set     => $alt_fset,
        });
};
ok($@, 'DataSet::reset_relational_attributes no -feature_set error');

eval { $clone_dset->reset_relational_attributes(
         {  
          -supporting_sets => $alt_ssets,
         })
};
ok($@, 'DataSet::reset_relational_attributes no -supporting_sets error');



#todo _set_Sets_and_types

done_testing();


1;

