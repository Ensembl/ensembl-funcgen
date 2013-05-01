#!usr/bin/env perl

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
#my $user = undef;

my $efgdba = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
    -user       => $user,
#    -pass       => ,
    -DNADB_USER => $user,
    #-DNADB_PORT => 3306,

    -species    => 'homo_sapiens',
    -dbname     => undef,
    -host       => undef,
    -DNADB_HOST => undef,
    -DNADB_NAME => 'homo_sapiens_core_71_37',
);

my $rsa = $efgdba->get_adaptor("resultset");


#Just grab a few result sets to work with
my ($result_set, $result_set_2) = 
  @{ $rsa->fetch_all_by_feature_class('dna_methylation', 
                                      {status => 'DISPLAYABLE'}) };
if(! (defined $result_set && defined $result_set_2)){
  throw('Failed to fetch 2 ResultSets to test, please update ResultSet.t');  
}


# START testing compare_to

my %diffs = %{$result_set->compare_to($result_set)};
my $diffs = '';
$diffs = "\n".Dumper(\%diffs) if %diffs;
ok(! %diffs, 'ResultSet::compare_to self defaults returns no diffs'.$diffs);

# redefine methods
%diffs = %{$result_set->compare_to($result_set, undef,
  [qw(name table_name feature_class get_all_states)] ,
  [qw(feature_type cell_type analysis get_support)]  
  )};
$diffs = "\n".Dumper(\%diffs) if %diffs;
ok(! %diffs, 'ResultSet::compare_to self redefined methods returns no diffs'.$diffs);

eval{#Invalid redefined scalar method
 $result_set->compare_to($result_set, undef, 
                   [qw(feature_type cell_type analysis get_support)] ); 
};
ok($@, 'Storable/ResultSet::compare_to scalar methods redefined (invalid)');


eval{#Invalid redefined object method
 $result_set->compare_to($result_set, undef, undef,
                   [qw(name table_name feature_class get_all_states)] ); 
};
ok($@, 'Storable/ResultSet::compare_to object methods redefined (invalid)');

# COMPLETED testing compare_to



## START testing reset_relational_attributes 


#clone ResultSet, so we can change some attrs
my $clone_rset   = bless({%{$result_set}}, ref($result_set));
my $alt_ftype    = $result_set->feature_type->adaptor->fetch_by_name('H4K4me3');
my $alt_analysis = $result_set->analysis->adaptor->fetch_by_logic_name('AFFY_UTR_ProbeAlign');
my $alt_ctype    = ($result_set_2->cell_type->name eq 'H1ESC') ?
                    $result_set_2->cell_type->adaptor->fetch_by_name('GM06990') :
                    $result_set_2->cell_type->adaptor->fetch_by_name('H1ESC');
my $alt_support  = $result_set_2->get_support;
my $orig_support = $result_set->get_support;#also lazyload now before we blow away the adaptor
#todo check support is different!

my %relational_params = 
  (
    -support      => $alt_support,
    -analysis     => $alt_analysis,
    -feature_type => $alt_ftype,
    -cell_type    => $alt_ctype,
  );

$clone_rset->reset_relational_attributes(\%relational_params, 'no_db_reset');


#eq comparisons of obj in scalar context compares mem refs
ok( ($clone_rset->analysis eq $alt_analysis),
   'ResultSet::reset_relational_attributes reset analysis');


#to test support ID differ here

ok( ($clone_rset->feature_type eq $alt_ftype),
   'ResultSet::reset_relational_attributes reset feature_type');

ok( ($clone_rset->cell_type eq $alt_ctype),
   'ResultSet::reset_relational_attributes reset ctype');
   
ok((defined $clone_rset->dbID && 
   defined $clone_rset->adaptor), 
   'ResultSet::reset_relational_attributes no_db_reset');
    
$clone_rset->reset_relational_attributes(\%relational_params);
ok(! (defined $clone_rset->dbID || 
      defined $clone_rset->adaptor), 
      'ResultSet::reset_relational_attributes with dbID/adaptor reset');


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



# COMPLETED testing reset_relational_attributes



# START testing add_support
#my @orig_support   = @{$orig_support};
#my @orig_support_2 = @{$result_set_2->get_support};
#push @orig_support, @orig_support_2;
my @new_support    = @{$clone_rset->add_support($orig_support)};#This appends rather than resets!

#let's assume if we have the same size, the we have been successful?
ok( ( (scalar(@new_support)) == 1 &&
      (scalar(@$orig_support))),
    'ResultSet::add_support returns expected size array');
    
ok( ( ($new_support[0] eq $orig_support->[0])),
    'ResultSet::add_support returns expected object');
    
    
    

#Duplicate support test has to be done on result_set_2 
#to avoid borking the rest of the tests
eval { $clone_rset->add_support($orig_support) };
ok($@, 'ResultSet::add_support caught duplicate support addition');  
   
   
  


done_testing();