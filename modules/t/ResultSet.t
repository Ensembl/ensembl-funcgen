#!usr/bin/env perl

use strict;
use warnings;

use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Test::More;
use Data::Dumper qw(Dumper);
use Bio::EnsEMBL::Utils::Exception qw( throw );

#throw('Test DB not yet implemented, you need to define a DBAdaptor and remove this throw manually');
my $user = '';

my $efgdba = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
    -user       => $user,
#    -pass       => ,
    -DNADB_USER => $user,
    #-DNADB_PORT => 3306,

    -species    => 'homo_sapiens',
    #-dbname     => 'nj1_test_homo_sapiens_funcgen_71_37',
    -host       => ,
    -DNADB_HOST => ,
    #-DNADB_NAME => 'homo_sapiens_core_71_37',
);



my $rsa = $efgdba->get_adaptor("resultset");

#Just grab a few result sets to work with
my ($result_set, $result_set_2) = 
  @{ $rsa->fetch_all_by_feature_class('dna_methylation', 
                                      {status => 'DISPLAYABLE'}) };
if(! (defined $result_set && defined $result_set_2)){
  throw('Failed to fetch 2 ResultSets to test, please update ResultSet.t');  
}

## START testing reset_relational_attributes and compare_to                                                 
#TODO Make sure $diffs only reports diffs for that test, 
#i.e. do we have any older diffs hanging around from previous tests

#Need lazy load via get_all_states before we blow away the adaptor
$result_set->get_all_states;
my $diffs = '';
my %diffs = %{$result_set->compare_to($result_set)};
$diffs = "\n".Dumper(\%diffs) if %diffs;
ok(! %diffs, 'ResultSet::compare_to self diffs'.$diffs);

#This is actually the same as above, just omiting the compare_stored_Storable
#checks on the nested objects/methods
%diffs = %{$result_set->compare_to($result_set, 'shallow')};
$diffs = "\n".Dumper(\%diffs) if %diffs;
ok(! %diffs, 'ResultSet::compare_to shallow compare_to self'.$diffs);

#clone ResultSet, so we can change some attrs
my $clone_rset = bless({%{$result_set}}, ref($result_set));
my %relational_params = 
  (
    -support      => $result_set->get_support,
    -analysis     => $result_set->analysis,
    -feature_type => $result_set->feature_type,
    -cell_type    => $result_set->cell_type,
  );

$clone_rset->reset_relational_attributes(\%relational_params, 'no_db_reset');
%diffs = %{$result_set->compare_to($clone_rset)};
$diffs = "\n".Dumper(\%diffs) if %diffs;
ok(! %diffs, 
   'ResultSet::compare_to clone after reset_relational_attributes with no change'
   .$diffs);
   
ok((defined $clone_rset->dbID && 
   defined $clone_rset->adaptor), 
   'ResultSet::reset_relational_attributes no db reset');
    
$clone_rset->reset_relational_attributes(\%relational_params);
ok(! (defined $clone_rset->dbID || 
      defined $clone_rset->adaptor), 
      'ResultSet::reset_relational_attributes with db reset');
  
%diffs = %{$result_set->compare_to($clone_rset)};
$diffs = "\n".Dumper(\%diffs) if %diffs;
ok(! %diffs, 
  'ResultSet::compare_to clone after reset_relational_attributes with db reset'
  .$diffs);

eval { $result_set->compare_to('NOT A ResultSet') };
ok($@, 'ResultSet::compare_to validate ResultSet arg');

eval { $clone_rset->reset_relational_attributes(
         {  
          -support      => $result_set->get_support,
          -analysis     => $result_set->analysis,
          -feature_type => $result_set->feature_type,
         })
};
ok($@, 'ResultSet::reset_relational_attributes no -cell_type error');

eval { $clone_rset->reset_relational_attributes(
         {  
          -support   => $result_set->get_support,
          -analysis  => $result_set->analysis,
          -cell_type => $result_set->cell_type,
         })
};
ok($@, 'ResultSet::reset_relational_attributes no -feature_type error');

eval { $clone_rset->reset_relational_attributes(
         {  
          -support   => $result_set->get_support,
          -cell_type => $result_set->cell_type,
          -feature_type => $result_set->feature_type,
         })
};
ok($@, 'ResultSet::reset_relational_attributes no -analysis error');

eval { $clone_rset->reset_relational_attributes(
         {  
          -cell_type => $result_set->cell_type,
          -feature_type => $result_set->feature_type,
          -analysis => $result_set->analysis,
         })
};
ok($@, 'ResultSet::reset_relational_attributes no -support error');

#Now test we identify different attributes
#Only need to change 1 of each of the Storable tested methods (string/objects)
#as those tests should be done in Storable.t not here
#we only want to test that we are catching the return value correctly

#Get the alternate objects for replacement
#These needs to be stored for compare_stored_Storables test
#my $alt_ctype = ($result_set->cell_type eq 'H1ESC') ? 
#  $result_set->cell_type->adaptor->fetch_by_name('GM06990') :
#  $result_set->cell_type->adaptor->fetch_by_name('H1ESC');

#We know it's already a DNA methylation feature typeof some sort
my $alt_ftype    = $result_set->feature_type->adaptor->fetch_by_name('H4K4me3');
#my $alt_analysis = $result_set->analysis->adaptor->fetch_by_logic_name('AFFY_UTR_ProbeAlign');

#todo reset all relation attributes at once and test compare_to for all in one go
#like FeatureSet.t

$clone_rset->reset_relational_attributes(
      {  
        -cell_type => $result_set->cell_type,
        -support   => $result_set->get_support,
        -feature_type => $alt_ftype,
        -analysis => $result_set->analysis,
      });

#Test passing diffs hash and catch ftype diff at same time
%diffs = %{$result_set->compare_to($clone_rset)};
ok( ( exists $diffs{'feature_type'} &&
     ($diffs{'feature_type'} =~ /^dbID mismatch/) 
    ), 
    'ResultSet::compare_to caught different feature_type '.
      '(compare_stored_Storables)');
  
# COMPLETED testing reset_relational_attributes  
$clone_rset->{name} = 'TEST_NAME';
%diffs = %{$result_set->compare_to($clone_rset)};
ok(exists $diffs{'name'}, 
  'ResultSet::compare_to caught different name '.
    '(compare_string_methods)');
$clone_rset->{name} = $result_set->name;

# START testing add_support and compare_to support
my @orig_support   = @{$clone_rset->get_support};
my @orig_support_2 = @{$result_set_2->get_support};
push @orig_support, @orig_support_2;
my @new_support    = @{$clone_rset->add_support(\@orig_support_2)};#This appends rather than resets!

#let's assume if we have the same size, the we have been successful?
ok( (scalar(@new_support) == scalar(@orig_support)),
    'ResultSet::add_support returns expected size array');

#Dupilcate support test has to be done on result_set_2 
#to avoid borking the rest of the tests
eval { $result_set_2->add_support(\@orig_support_2) };
ok($@, 'ResultSet::add_support caught duplicate support addition');  
   
#TODO add more test on add_support?
#should really test we have exactly the same support objects

#This now fails as we don't do the size test in shallow mode. 
%diffs = %{$result_set->compare_to($clone_rset)};
ok( (exists $diffs{'get_support'} &&
    ($diffs{'get_support'} =~ /Return size mismatch/) ), 
   'ResultSet::compare_to support size');
    
$clone_rset->reset_relational_attributes(
      {  
        -cell_type => $result_set->cell_type,
        -support   => $result_set_2->get_support,
        -feature_type => $result_set->feature_type,
        -analysis => $result_set->analysis,
      });    
#This will return no diffs even tho the dbID differ as we are only doing a shallow compare   
%diffs = %{$result_set->compare_to($clone_rset, 'shallow')};
$diffs = "\n".Dumper(\%diffs) if %diffs;
ok(! %diffs,'ResultSet::compare_to shallow overlook support diffs'.$diffs);

%diffs = %{$result_set->compare_to($clone_rset)};
ok( (exists $diffs{'get_support'} &&
     ($diffs{'get_support'} =~ /dbID mismatch/) ),
  'ResultSet::compare_to get_support non shallow dbID mismatch');

#Finished testing compare_to

done_testing();
