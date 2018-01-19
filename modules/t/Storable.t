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

use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Test::More;
use Data::Dumper qw(Dumper);
use Bio::EnsEMBL::Utils::Exception qw( throw );

# This test uses the following comment conventions
# START method_name testing
# COMPLETED method_name testing

throw('Test DB not yet implemented, you need to define a DBAdaptor and remove this throw manually');
my $user = undef;

my $efgdba = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
    -user       => $user,
#    -pass       => ,
    -DNADB_USER => $user,
    #-DNADB_PORT => 3306,

    -species    => 'homo_sapiens',
    -host       => undef,
    -dbname     => undef,
    -DNADB_HOST => undef,
    -DNADB_NAME => 'homo_sapiens_core_71_37',
);


my $error;

#Just grab a few result sets to work with
#DISPLAYABLE ensure we should have an input_set defined
my $fsa = $efgdba->get_adaptor('featureset');
my ($fset, $fset_2) =  @{ $fsa->fetch_all_by_feature_class('annotated', 'DISPLAYABLE')};
if(! (defined $fset && defined $fset_2)){
  throw('Failed to fetch 2 FeatureSets to test, please update FeatureSet.t');  
}

my $ftype = $fset->feature_type;


#How do we split these up wrt dataclass and storable tests
#can't easily do all test in Storable as this won't have the
#required methods to test


## START testing compare_to compare_object_methods compare_scalar_methods
#_compare_method_return_types (implicit through above)
# using FeatureSet::reset relational_relational_attributes  

#### Test compare_to mandatory params

eval{#No storable arg  
  $fset->compare_to();
};
ok($@, 'Storable::compare_to Storable not defined');

eval{#Different class 
 $fset->compare_to($ftype); 
};
ok($@, 'Storable::compare_to same class');


eval{#Invalid redefined scalar method
 $fset->compare_to($fset, undef, 
                   [qw(feature_type cell_type analysis get_InputSet)] ); 
};
ok($@, 'Storable/FeatureSet::compare_to scalar methods redefined (invalid)');


eval{#Invalid redefined object method
 $fset->compare_to($fset, undef, undef,
                   [qw(name description display_label feature_class get_all_states)] ); 
};
ok($@, 'Storable/FeatureSet::compare_to object methods redefined (invalid)');

#TODO add tests for array return types (and arrayref for objects)
#Need to change this to ResultSet tests so we can check size mismatch here too for objs?
#what will we actually test in the data classes?  
#Re-iterate the simply test of calling with method args then swapping the method defs

                                              
#Make sure we lazy load via get_all_states before we blow away the adaptor
#$fset->get_all_states; #This will be called in the compare_to
my $diffs = '';
my %diffs = %{$fset->compare_to($fset)};
$diffs = "\n".Dumper(\%diffs) if %diffs;
ok(! %diffs, 'Storable::compare_to self diffs'.$diffs);

#re-define the methods explicitly
%diffs = %{$fset->compare_to($fset, undef,
  [qw(name description display_label feature_class get_all_states)] ,
  [qw(feature_type cell_type analysis get_InputSet)] 
  )};
$diffs = "\n".Dumper(\%diffs) if %diffs;
ok(! %diffs, 'Storable::compare_to compare_to self, re-defined methods'.$diffs);


#This is actually the same as above, 
#just omiting the nested object dbID/adaptor
%diffs = %{$fset->compare_to($fset, 1)};
$diffs = "\n".Dumper(\%diffs) if %diffs;
ok(! %diffs, 'Storable::compare_to shallow compare_to self'.$diffs);

#clone FeatureSet, so we can change some attrs
my $clone_fset = bless({%{$fset}}, ref($fset));
 
#Make sure we don't catch the dbID/adaptor change
undef $clone_fset->{dbID};  
undef $clone_fset->{adaptor};
%diffs = %{$fset->compare_to($clone_fset)};
$diffs = "\n".Dumper(\%diffs) if %diffs;
ok(! %diffs, 
  'Storable::compare_to clone without dbID/adaptor'
  .$diffs);


#Test invalid compare_to arg
eval { $fset->compare_to('NOT A FeatureSet') };
ok($@, 'Storable::compare_to validate FeatureSet arg');


my $alt_iset = $fset_2->get_InputSet;
if( (! $alt_iset) || 
    ($alt_iset->dbID == $fset->get_InputSet->dbID) ){
  throw('Failed to fetch alternative InputSet for reset_relational_attributes.'.
          "\nPlease amend FeatureSet.t"); 
}

#Only need to change one here as we are actually testing
#that compare_to picks up an object_method change

$clone_fset->reset_relational_attributes(
      {  
        -cell_type    => $fset->cell_type,
        -input_set    => $alt_iset,
        -feature_type => $fset->feature_type,
        -analysis     => $fset->analysis,
      });

#Test passing diffs hash and catch ftype diff at same time
%diffs = %{$fset->compare_to($clone_fset)};

ok((exists $diffs{'get_InputSet'} &&
     ($diffs{'get_InputSet'} =~ /^dbID mismatch/) ),
  'FeatureSet::reset_relational_attributes and Storable::compare_to caught different input_set '.
    '(object_methods)');    
 
#Don't need to do other string methods here, as we are testing
#that the return of compare_string method is caught properly, 
#not the individual string methods
$clone_fset->{name} = 'TEST_NAME';
%diffs = %{$fset->compare_to($clone_fset)};
ok(exists $diffs{'name'}, 
  'Storable::compare_to/compare_scalar_methods caught different name '.
    '(scalar_methods) in passed hashref');
$clone_fset->{name} = $fset->name;

#This should return no diffs even though we have a alt_iset
%diffs = %{$fset->compare_to($clone_fset, 1)};
$diffs = "\n".Dumper(\%diffs) if %diffs;
ok(! %diffs,'Storable::compare_to shallow overlooks object_methods'.$diffs);


#START testing has_status add_status
my $status = 'RESOLVED';#fset should never had resolved
$clone_fset->{states} = [];
ok( ! $clone_fset->has_status($status),
    'Storable::has_status returns false'); 

eval{#Mandatory param
  $clone_fset->add_status(); 
};
ok($@,
   'Storable::add_status mandatory parameter');

$clone_fset->add_status($status);
ok( $clone_fset->has_status($status),
    'Storable::add_status and has_status return true');
 
#todo do we need to test fetching states from DB here? 
#COMPLETED testing add_status has_status

#compare_method_return_types arrayref
eval { 
  #This may throw an error if it doesn't work properly
  #Have to specify empty listref to override default object methods
  %diffs = %{$fset->compare_to($clone_fset, undef, 
                             [qw(get_all_states)], [])};
};

$error = $@;
ok(( (! $error) && 
     (scalar(keys %diffs) == 1) &&
     (exists $diffs{get_all_states}) ),
  'Storable::_compare_method_return_types arrayref '.$error);
 
#Storable::compare_scalar_methods caught non-scalar return types
#Need to reset this in both as method will assume same return types
#and only check one object($fset)
my $tmp_ftype = $fset->feature_type; #for restore later?
my $test_hash = {};
$clone_fset->{states} = [$test_hash];
$fset->{states}       = [$test_hash];
 
eval { 
  #Have to specify empty listref to override default object methods
  %diffs = %{$fset->compare_to($clone_fset, undef, 
                             [qw(get_all_states)], [])};
};


ok( $@,
  'Storable::compare_scalar_methods caught non-scalar return types');
 
 
 
 
#compare_method_return_types caught returned ref is not ARRAY
$clone_fset->{feature_type} = $test_hash;
$fset->{feature_type}       = $test_hash;

eval{
  #Have to specify empty listref to override default scalar methods
  %diffs = %{$fset->compare_to($clone_fset, undef, [],
                               [qw(feature_type)])};
};
ok( $@,
  'Storable::_compare_method_return_types caught returned ref is not ARRAY');
 
 
#Storable::compare_storable_methods caught return type is not a Storable
$clone_fset->{feature_type} = [$test_hash, $test_hash];
$fset->{feature_type}       = [$test_hash, $test_hash];
eval{
  #Have to specify empty listref to override default scalar methods
  %diffs = %{$fset->compare_to($clone_fset, undef, [],
                               [qw(feature_type)])};
};
ok( $@,
  'Storable::_compare_method_return_types caught storable does not return Storable(s)');
 

#Storable::compare_storable_methods caught Storables namespace mismatch
$clone_fset->{feature_type} = $fset;
$fset->{feature_type}       = $tmp_ftype;

eval{
  #Have to specify empty listref to override default scalar methods
  %diffs = %{$fset->compare_to($clone_fset, undef, [],
                               [qw(feature_type)])};
};

$error = $@;

ok( ((!$error) && 
     (scalar(keys %diffs) == 1) &&
     (exists $diffs{feature_type}) && 
     ($diffs{feature_type} =~ /Namespace mismatch/)),
  'Storable::compare_storable_methods caught Storables namespace mismatch '.$error);
 



#we can't fake the array return type as this would requires modification of the method
#do we have any available method that returns and array?

#todo add test for undef returned from both scalar/storable
#todo add test for objects from different DBs
#This is quite hard to do as we will need to have two different objects with the same method name which 
#return different classes. Fake this up!

# COMPLETED testing compare_to compare_object_methods compare_scalar_methods
#_compare_method_return_types (implicit through above)



done_testing();



