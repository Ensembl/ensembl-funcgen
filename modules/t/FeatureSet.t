#!usr/bin/env perl

use strict;
use warnings;

use Test::More;
use Data::Dumper                   qw( Dumper );
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Utils::Scalar    qw (check_ref );
#use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;


# This test uses the following comment conventions
# START method_name testing
# COMPLETED method_name testing

throw('Test DB not yet implemented, you need to define a DBAdaptor and remove this throw manually');

my $efgdba = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
    -user       => 'XXX',
    #-DNADB_USER => $dnadb_user,
    -species    => 'homo_sapiens',
    -dbname     => 'homo_sapiens_funcgen_78_38',
    -host       => 'XXX',
    #-DNADB_HOST => $dnadb_host,
    #-DNADB_NAME => $dnadb_name,
);
$efgdba->dbc->db_handle;#Test DB

my $fsa = $efgdba->get_adaptor("featureset");

#Just grab a few result sets to work with
#DISPLAYABLE ensure we should have an input_set defined
my ($fset, $fset_2) =  @{ $fsa->fetch_all_by_feature_class('annotated', 'DISPLAYABLE')};
if(! (defined $fset && defined $fset_2)){
  throw('Failed to fetch 2 FeatureSets to test, please update FeatureSet.t');  
}

# START testing compare_to

# TODO
# Need to eval and test these compare_to calls, as they may fail, abroting the test early.
# This can happen if the API changes and the methods passed here are no longer valid

my $diffs = '';
my %diffs = %{$fset->compare_to($fset)};
$diffs = "\n".Dumper(\%diffs) if %diffs;
ok(! %diffs, 'FeatureSet::compare_to self default nodiffs'.$diffs);

# redefine methods
%diffs = %{$fset->compare_to($fset, undef,
  [qw(name description display_label feature_class get_all_states)] ,
  [qw(feature_type cell_type analysis)]  
  )};
$diffs = "\n".Dumper(\%diffs) if %diffs;
ok(! %diffs, 'FeatureSet::compare_to self redefined methods no diffs'.$diffs);

eval{#Invalid redefined scalar method
 $fset->compare_to($fset, undef, 
                   [qw(feature_type cell_type analysis)] ); 
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

my $alt_exp = $fset_2->experiment;

if( (! $alt_exp) || 
    ($alt_exp->dbID == $fset->experiment->dbID) ){
  throw('Failed to fetch alternative Experiment for reset_relational_attributes.'.
          "\nPlease amend FeatureSet.t or test DB"); 
}


SKIP: {

  if(! ($alt_ctype && $alt_ftype && $alt_anal)){
    skip("Fkipping reset_relational_attributes with alterates test as failed to acquire alternative relational attributes:\n\t".
           "Analysis:\t${alt_anal}\n\tCellType:\t${alt_ctype}\n\tFeatureType\t${alt_ftype}".
           "\nPlease amend FeatureSet.t or test DB", 5); 
  }

  #print "Using alternate attributes:\t".join("\t", ($alt_ctype->name, $alt_ftype->name, $alt_anal->logic_name))."\n";
  my %relational_params = 
   (-experiment   => $alt_exp,
    -analysis     => $alt_anal,
    -feature_type => $alt_ftype,
    -cell_type    => $alt_ctype);

  #Need to eval this in case of failure
  $clone_fset->reset_relational_attributes(\%relational_params, 'no_db_reset');

  #eq comparisons of obj in scalar context compares mem refs
  ok(($clone_fset->analysis eq $alt_anal),
     'FeatureSet::reset_relational_attributes reset analysis');

  ok(($clone_fset->experiment eq $alt_exp),
     'FeatureSet::reset_relational_attributes reset experiment');

  ok(($clone_fset->feature_type eq $alt_ftype),
     'FeatureSet::reset_relational_attributes reset feature_type');

  ok(($clone_fset->cell_type eq $alt_ctype),
     'FeatureSet::reset_relational_attributes reset ctype');

  ok( ( (defined $clone_fset) && (defined $clone_fset->adaptor) ),
     'FeatureSet::reset_relational_attributes no_db_reset');
   
  $clone_fset->reset_relational_attributes(\%relational_params);
  ok( ((! defined $clone_fset->dbID) && (! defined $clone_fset->adaptor) ),
     'FeatureSet::reset_relational_attributes with dbID/adaptor reset');
}

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
          -experiment => $fset->experiment,
          -analysis   => $fset->analysis,
          -cell_type  => $fset->cell_type,
         })
};
ok($@, 'FeatureSet::reset_relational_attributes no -feature_type error');

eval { $clone_fset->reset_relational_attributes(
         {  
          -experiment   => $fset->experiment,
          -cell_type    => $fset->cell_type,
          -feature_type => $fset->feature_type,
         })
};
ok($@, 'FeatureSet::reset_relational_attributes no -analysis error');

eval { $clone_fset->reset_relational_attributes(
         {  
          -cell_type    => $fset->cell_type,
          -feature_type => $fset->feature_type,
          -analysis     => $fset->analysis,
         })
};
ok($@, 'FeatureSet::reset_relational_attributes no -experiment error');

#Test optional reset_relational_attributes params
#conditional on them not being present in the object i.e. undef them first
undef $clone_fset->{experiment};
undef $clone_fset->{experiment_id};
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
ok(! $@, 'FeatureSet::reset_relational_attributes optional -experiment');
          
undef $clone_fset->{cell_type};
eval { $clone_fset->reset_relational_attributes(
         {  
          -feature_type => $fset->feature_type,
          -analysis     => $fset->analysis,
         })
};
ok(! $@, 'FeatureSet::reset_relational_attributes optional -cell_type');

# COMPLETED testing reset_relational_attributes  



#Mouse specific test 

my $mdb = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
           -user       => 'XXX',
           -species    => 'mus_musculus',
           -dbname     => 'mus_musculus_funcgen_78_38',
           -host       => 'XXX',
);

my $mirna = 'mmu-miR-16-5p';#MIMAT0000527
my $ftype = $mdb->get_FeatureTypeAdaptor->fetch_by_name($mirna);
my $fset_name = 'TarBase miRNA';
$fset  = $mdb->get_FeatureSetAdaptor->fetch_by_name($fset_name);

SKIP: {
  if(! ($ftype && $fset)){
    skip('Skipping get_Features_by_FeatureType test. '.
     "Failed to fetch $fset_name FeatureSet and/or $mirna FeatureType", 3);
  }

  my $mirna_feats = $fset->get_Features_by_FeatureType($ftype);

  my $is_aref = check_ref($mirna_feats, 'ARRAY');

  ok($is_aref, 'get_Features_by_FeatureType returns Arrayref');


  SKIP: {
    if(! $is_aref){
      skip('Skipping get_Features_by_FeatureType count test. Failed to return Arrayref', 2);
    }


    #should be 1019 in mus_musculus_funcgen_78_38
    my $true_ft_cnt = 1019;
    my $feat_cnt = scalar(@$mirna_feats);
    ok($feat_cnt == $true_ft_cnt, "Got $true_ft_cnt miRNATargetFeatures via get_Features_by_FeatureType:\t$feat_cnt");

    #Test all have correct FeatureType.

    my @mismatches = grep {! /$mirna/ } (map {$_->feature_type->name} @$mirna_feats);
    ok(! @mismatches, "Found no mismatches FeatureTypes from get_Features_by_FeatureType:\t".scalar(@mismatches));  
  }


}


done_testing();



