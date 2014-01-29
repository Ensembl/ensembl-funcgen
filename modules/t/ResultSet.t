#!usr/bin/env perl

use strict;
use warnings;

use Test::More                     qw( no_plan );
use Data::Dumper                   qw( Dumper );
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Utils::Scalar    qw( check_ref );
use Bio::EnsEMBL::Funcgen::ResultSet;
#use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

#obtain Adaptors for dnabb and funcgen databases
#my $multi  = Bio::EnsEMBL::Test::MultiTestDB->new();
#my $efgdba = $multi->get_DBAdaptor("funcgen");

# This test uses the following comment conventions
# START method_name testing
# COMPLETED method_name testing

throw('MultiTestDB not yet implemented, you need to temporarily define a DBAdaptor and remove this throw manually');
my $user       = undef;
my $dnadb_user = $user;
my $host       = undef;
my $dbname     = undef;
my $dnadb_host = undef;
my $dnadb_name = undef;

my $efgdba = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
    -user       => $user,
    -DNADB_USER => $dnadb_user,
    #-species    => 'homo_sapiens',
    -dbname     => $dbname,
    -host       => $host,
    -DNADB_HOST => $dnadb_host,
    -DNADB_NAME => $dnadb_name,
);
$efgdba->dbc->db_handle;#Test DB

my $rsa      = $efgdba->get_adaptor("resultset");
my $slice_a  = $efgdba->dnadb->get_SliceAdaptor;
my $exp_name = 'H1ESC_Tcf12_ENCODE_Hudsonalpha';
my $exp      = $efgdba->get_adaptor("experiment")->fetch_by_name($exp_name);


    
SKIP:{
  if(! $exp){
    skip "Could not fetch test Experiment:\t$exp_name\nThis must have been removed from the DB, please choose another Experiment or fix the test DB";
  }
    
  my @rsets;

  eval { @rsets = @{$rsa->fetch_all_by_Experiment($exp)}; };
  ok(! $@, "ResultSetAdaptor::fetch_all_by_Experiment executes without error:\t$@");

  ok(scalar(@rsets) == 1,
   'ResultSetAdaptor::fetch_all_by_Experiment returns expected number of ResultSets');

  #We now append the analysis logic name, so let's strip that off for our test
  #This is heavily dependant on our usage of ResultSets and the naming scheme
  my $lname          = $rsets[0]->analysis->logic_name;
  my $rset_name      = $rsets[0]->name;
  (my $rset_exp_name = $rset_name) =~ s/_${lname}//;

  ok($rset_exp_name eq $exp_name,
   "ResultSetAdaptor::fetch_all_by_Experiment returns expected ResultSet: Experiment $exp_name is link to ResultSet ".$rset_name);
   
   
  #Now let's get and check the file path for each window size
  #This really spans ResultSet and ResultFeatureAdaptor
  
  #todo We need to check for dbfile.data_root firt and skip
  
  #This is actually tested in the CheckResultSetDBFileLinks HC
  #but here we are checking the API, not the data
  #although we are using the fact that the data file exists to validate the API method

  my $slice = $slice_a->fetch_by_region('chromosome', 1, 100000, 1100000);

  
  foreach my $wsize(@{$efgdba->get_ResultFeatureAdaptor->window_sizes}){#8 wsize tests here
    my $col_file = $rsets[0]->get_dbfile_path_by_window_size($wsize);
    ok(-f $col_file, "Collection (window_size=$wsize) file for $rset_name exists:\n".
      "\t$col_file");    
    
    #Now let's test the get_ResultFeatures_by_Slice wrapper for each window?
    my $rfs = $rsets[0]->get_ResultFeatures_by_Slice($slice, undef, $wsize);
   
    #Let's assume that the method wrapped is functioning properly
    #so we don't have to recreate the tests for that
    #basically we want to test that the args are bing passed on in the correct way
    #would require recreating nested skips for the return types before we can validate the wsize
    my $rtype_valid = check_ref($rfs, 'ARRAY');
    my $rf_valid    = check_ref($rfs->[0], 'Bio::EnsEMBL::Funcgen::Collection::ResultFeature');
    my $rf_wsize;
    eval { $rf_wsize = $rfs->[0]->window_size }; 
    
    ok($rtype_valid && (scalar(@$rfs) == 1) && 
       $rf_valid && ($rfs->[0]->window_size == $wsize),
       'Got expected single ResultFeature Collection '.
       "(window_size=$wsize) for $rset_name:\t@$rfs (window_size=$rf_wsize)"); 
  }
}


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
