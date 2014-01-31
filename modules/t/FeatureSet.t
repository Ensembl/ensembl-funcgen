#!usr/bin/env perl

use strict;
use warnings;

use Test::More;
use Data::Dumper                   qw( Dumper );
use Bio::EnsEMBL::Utils::Exception qw( throw );
#use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;


# This test uses the following comment conventions
# START method_name testing
# COMPLETED method_name testing

throw('Test DB not yet implemented, you need to define a DBAdaptor and remove this throw manually');

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

my $rsa = $efgdba->get_adaptor("featureset");

#Just grab a few result sets to work with
#DISPLAYABLE ensure we should have an input_set defined
my ($fset, $fset_2) =  @{ $rsa->fetch_all_by_feature_class('annotated', 'DISPLAYABLE')};
if(! (defined $fset && defined $fset_2)){
  throw('Failed to fetch 2 FeatureSets to test, please update FeatureSet.t');  
}

# START testing compare_to

my $diffs = '';
my %diffs = %{$fset->compare_to($fset)};
$diffs = "\n".Dumper(\%diffs) if %diffs;
ok(! %diffs, 'FeatureSet::compare_to self default nodiffs'.$diffs);

# redefine methods
%diffs = %{$fset->compare_to($fset, undef,
  [qw(name description display_label feature_class get_all_states)] ,
  [qw(feature_type cell_type analysis get_InputSet)]  
  )};
$diffs = "\n".Dumper(\%diffs) if %diffs;
ok(! %diffs, 'FeatureSet::compare_to self redefined methods no diffs'.$diffs);

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


# COMPLETED testing compare_to



## START testing reset_relational_attributes                                             

#clone ResultSet, so we can change some attrs
my $clone_fset = bless({%{$fset}}, ref($fset));

#Get the alternate objects for replacement
my $alt_ftype = $fset->feature_type->adaptor->fetch_by_name('RegulatoryFeature');
my $alt_ctype = ($fset->cell_type->name eq 'H1ESC') ?
  $fset->cell_type->adaptor->fetch_by_name('GM06990') :
  $fset->cell_type->adaptor->fetch_by_name('H1ESC');
my $alt_anal =  ($fset->analysis->logic_name eq 'ccat_histone') ?
  $fset->analysis->adaptor->fetch_by_logic_name('AFFY_UTR_ProbeAlign') :
  $fset->analysis->adaptor->fetch_by_logic_name('ccat_histone'); 
  
if(! ($alt_ctype && $alt_ftype && $alt_anal)){
  throw("Failed to acquire alternative relational attributes:\n\t".
          "Analysis:\t${alt_anal}\n\tCellType:\t$alt_ctype\n\tFeatureType\t${alt_ftype}".
          "\nPlease amend FeatureSet.t"); 
}

my $alt_iset = $fset_2->get_InputSet;

if( (! $alt_iset) || 
    ($alt_iset->dbID == $fset->get_InputSet->dbID) ){
  throw('Failed to fetch alternative InputSet for reset_relational_attributes.'.
          "\nPlease amend FeatureSet.t"); 
}

my %relational_params = 
  (
   -input_set    => $alt_iset,
   -analysis     => $alt_anal,
   -feature_type => $alt_ftype,
   -cell_type    => $alt_ctype,
  );

$clone_fset->reset_relational_attributes(\%relational_params, 'no_db_reset');

#eq comparisons of obj in scalar context compares mem refs
ok( ($clone_fset->analysis eq $alt_anal),
   'FeatureSet::reset_relational_attributes reset analysis');

ok( ($clone_fset->get_InputSet eq $alt_iset),
   'FeatureSet::reset_relational_attributes reset input_set');

ok( ($clone_fset->feature_type eq $alt_ftype),
   'FeatureSet::reset_relational_attributes reset feature_type');

ok( ($clone_fset->cell_type eq $alt_ctype),
   'FeatureSet::reset_relational_attributes reset ctype');

ok( ( (defined $clone_fset) && (defined $clone_fset->adaptor) ),
   'FeatureSet::reset_relational_attributes no_db_reset');
   
$clone_fset->reset_relational_attributes(\%relational_params);
ok( ((! defined $clone_fset->dbID) && (! defined $clone_fset->adaptor) ),
   'FeatureSet::reset_relational_attributes with dbID/adaptor reset');

#Test reset_relational_attributes mandatory params
eval { $clone_fset->reset_relational_attributes(
         {
          -input_set    => $fset->input_set,   
          -analysis     => $fset->analysis,
          -feature_type => $fset->feature_type,
         })
};
ok($@, 'FeatureSet::reset_relational_attributes no -cell_type error');

eval { $clone_fset->reset_relational_attributes(
         {  
          -input_set => $fset->get_InputSet,
          -analysis  => $fset->analysis,
          -cell_type => $fset->cell_type,
         })
};
ok($@, 'FeatureSet::reset_relational_attributes no -feature_type error');

eval { $clone_fset->reset_relational_attributes(
         {  
          -input_set    => $fset->get_InputSet,
          -cell_type    => $fset->cell_type,
          -feature_type => $fset->feature_type,
         })
};
ok($@, 'FeatureSet::reset_relational_attributes no -analysis error');


#we have already tested that we have an input_set above
eval { $clone_fset->reset_relational_attributes(
         {  
          -cell_type    => $fset->cell_type,
          -feature_type => $fset->feature_type,
          -analysis     => $fset->analysis,
         })
};
ok($@, 'FeatureSet::reset_relational_attributes no -input_set error');

#Test optional reset_relational_attributes params
#conditional on them not being present in the object i.e. undef them first
undef $clone_fset->{input_set};
undef $clone_fset->{input_set_id};
#problem here is that we have removed the adaptor
#and get_InputSet tries to use the adaptor to fetch the input set if it not defined

eval { 
  $clone_fset->reset_relational_attributes(
         {  
          -feature_type => $fset->feature_type,
          -analysis     => $fset->analysis,
          -cell_type    => $fset->cell_type,
         });
};
ok(! $@, 'FeatureSet::reset_relational_attributes optional -input_set');
          
undef $clone_fset->{cell_type};
eval { $clone_fset->reset_relational_attributes(
         {  
          -feature_type => $fset->feature_type,
          -analysis     => $fset->analysis,
         })
};
ok(! $@, 'FeatureSet::reset_relational_attributes optional -cell_type');


# COMPLETED testing reset_relational_attributes  

done_testing();



