#!usr/bin/env perl

use strict;
use warnings;

use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Test::More;
use Data::Dumper qw(Dumper);
use Bio::EnsEMBL::Utils::Exception qw( throw );

throw('Test DB not yet implemented, you need to define a DBAdaptor and remove this throw manually');
my $user = '';

my $efgdba = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
    -user       => $user,
#    -pass       => ,
    -DNADB_USER => $user,
    #-DNADB_PORT => 3306,

    -species    => 'homo_sapiens',
    #-dbname     => 'nj1_test_homo_sapiens_funcgen_71_37',
    #-host       => 
    #-DNADB_HOST => 
    #-DNADB_NAME => 'homo_sapiens_core_71_37',
);



my $rsa = $efgdba->get_adaptor("featureset");

#Just grab a few result sets to work with
#DISPLAYABLE ensure we should have an input_set defined
my ($fset, $fset_2) =  @{ $rsa->fetch_all_by_feature_class('annotated', 'DISPLAYABLE')};
if(! (defined $fset && defined $fset_2)){
  throw('Failed to fetch 2 FeatureSets to test, please update FeatureSet.t');  
}

## START testing reset_relational_attributes and compare_to                                                 
#TODO Make sure $diffs only reports diffs for that test, 
#i.e. do we have any older diffs hanging around from previous tests

#Need lazy load via get_all_states before we blow away the adaptor
$fset->get_all_states;
my $diffs = '';
my %diffs = %{$fset->compare_to($fset)};
$diffs = "\n".Dumper(\%diffs) if %diffs;
ok(! %diffs, 'FeatureSet::compare_to self diffs'.$diffs);

#This is actually the same as above, just omiting the compare_stored_Storable
#checks on the nested objects/methods
%diffs = %{$fset->compare_to($fset, undef, 'shallow')};
$diffs = "\n".Dumper(\%diffs) if %diffs;
ok(! %diffs, 'FeatureSet::compare_to shallow compare_to self'.$diffs);

#clone ResultSet, so we can change some attrs
my $clone_fset = bless({%{$fset}}, ref($fset));
my %relational_params = 
  (
   -input_set    => $fset->get_InputSet,
   -analysis     => $fset->analysis,
   -feature_type => $fset->feature_type,
   -cell_type    => $fset->cell_type,
  );

#Test no diffs first
$clone_fset->reset_relational_attributes(\%relational_params, 'no_db_reset');
%diffs = %{$fset->compare_to($clone_fset)};
$diffs = "\n".Dumper(\%diffs) if %diffs;
ok(! %diffs, 
   'FeatureSet::compare_to clone after reset_relational_attributes with no change'
   .$diffs);
   
ok((defined $clone_fset->dbID && 
   defined $clone_fset->adaptor), 
   'FeatureSet::reset_relational_attributes no db reset');
    
$clone_fset->reset_relational_attributes(\%relational_params);
ok(! (defined $clone_fset->dbID || 
      defined $clone_fset->adaptor), 
      'FeatureSet::reset_relational_attributes with db reset');
  
%diffs = %{$fset->compare_to($clone_fset)};
$diffs = "\n".Dumper(\%diffs) if %diffs;
ok(! %diffs, 
  'FeatureSet::compare_to clone after reset_relational_attributes with db reset'
  .$diffs);


#Test invalid compare_to arg
eval { $fset->compare_to('NOT A FeatureSet') };
ok($@, 'FeatureSet::compare_to validate FeatureSet arg');

#This should never be true for annotated FeatureSets
if(! $fset->cell_type){
  throw('Attempting to test a FeatureSet with no CellType defined, please amend FeatureSet.t'); 
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


#This *should* never be true for DISPLAYABLE annotated FeatureSets
if(! $fset->get_InputSet){
  throw('Attempting to test a FeatureSet with no InputSet defined, please amend FeatureSet.t'); 
}
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


#Get the alternate object for replacement
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


#Test all reset_relational_attributes together with compare_to

#This is not catching an undef $alt_iset?

$clone_fset->reset_relational_attributes(
      {  
        -cell_type    => $alt_ctype,
        -input_set    => $alt_iset,
        -feature_type => $alt_ftype,
        -analysis     => $alt_anal,
      });

#Test passing diffs hash and catch ftype diff at same time
%diffs = %{$fset->compare_to($clone_fset)};
ok( (exists $diffs{'feature_type'} &&
     ($diffs{'feature_type'} =~ /^dbID mismatch/) ),
  'FeatureSet::reset_relational_attributes reset and compare_to caught different feature_type '.
    '(compare_stored_Storables)');
    
ok((exists $diffs{'cell_type'} &&
     ($diffs{'cell_type'} =~ /^dbID mismatch/) ),
  'FeatureSet::reset_relational_attributes reset and compare_to caught different cell_type '.
    '(compare_stored_Storables)');    
    
ok((exists $diffs{'analysis'} &&
     ($diffs{'analysis'} =~ /^dbID mismatch/) ),
  'FeatureSet::reset_relational_attributes reset and compare_to caught different analysis '.
    '(compare_stored_Storables)');    
    
ok((exists $diffs{'get_InputSet'} &&
     ($diffs{'get_InputSet'} =~ /^dbID mismatch/) ),
  'FeatureSet::reset_relational_attributes reset and compare_to caught different input_set '.
    '(compare_stored_Storables)');    
    
# COMPLETED testing reset_relational_attributes  

#Don't need to do other string methods here, as we are testing
#that the return of compare_string method is caguth properly, 
#not the individual string methods
$clone_fset->{name} = 'TEST_NAME';
%diffs = %{$fset->compare_to($clone_fset, \%diffs)};
ok(exists $diffs{'name'}, 
  'FeatureSet::compare_to caught different name '.
    '(compare_string_methods) in passed hashref');
$clone_fset->{name} = $fset->name;



#This should return no diffs 
%diffs = %{$fset->compare_to($clone_fset, 'shallow')};
$diffs = "\n".Dumper(\%diffs) if %diffs;
ok(! %diffs,'FeatureSet::compare_to shallow overlooks object methods'.$diffs);

#Finished testing compare_to

done_testing();



