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

my $rsa       = $efgdba->get_adaptor("resultset");
my $slice_a   = $efgdba->dnadb->get_SliceAdaptor;
my $rset_name = 'H1ESC_Tcf12_ENCODE_Hudsonalpha_bwa_samse';
my $rset      = $rsa->fetch_by_name($rset_name);
my $rf_a      = $efgdba->get_ResultFeatureAdaptor;
my $slice     = $slice_a->fetch_by_region('chromosome', 1, 100001, 1100000);
my $num_window_sizes = 8;


#All these SKIP blocks are to prevent exiting when trying to call a method on the return
#value of a previously failed test i.e. we don't want to call a method on an unblessed 
#reference which will cause the test to exit prematurely.


#todo test the negative cases here i.e. omit slice and rsets args
    
my ($rtype_valid, $rfs);
SKIP:{
  if(! $rset){
    skip "Could not fetch test ResultSet:\t$rset_name\nThis must have been removed from the DB.\n".
      'Please choose another Result or fix the test DB';
  }
   
  #Does this first so we don't have redundant tests in the loop
  $rfs = $rf_a->fetch_all_by_Slice_ResultSets($slice, [$rset]);   
  $rtype_valid = check_ref($rfs, 'HASH');    
  ok($rtype_valid && (scalar(keys %$rfs) == 1),
    'fetch_all_by_Slice_ResultSets returned a HASHREF containing a single key value pair:'.
    "\t$rfs (keys are:\t".join(' ', keys %$rfs).')');
}


my ($vtype_valid, $rfcs);
SKIP:{  
  skip "fetch_all_by_Slice_ResultSets return type not valid" if ! $rtype_valid;
  
  ($rfcs) = values %$rfs;
  #arbitrarily take the first if the above has failed to return a single pair
  $vtype_valid = check_ref($rfcs, 'ARRAY');
  ok($vtype_valid && (scalar(@$rfcs) == 1),
     'fetch_all_by_Slice_ResultSets returned a HASHREF value which is an ARRAYREF '.
     "with a single element:\t$rfcs");
}


my ($rf_ref_valid, $rfc);
SKIP:{  
  skip "fetch_all_by_Slice_ResultSets returned HASHREF values not valid" if ! $vtype_valid;
 
  #Arbitrarily take the first if the above failed to return a single array element
  ($rfc) = $rfcs->[0];
  $rf_ref_valid = check_ref($rfc, 'Bio::EnsEMBL::Funcgen::Collection::ResultFeature');
  ok($rf_ref_valid, 'fetch_all_by_Slice_ResultSets returned a HASHREF containing a '.
                      " ResultFeature Collection value:\t".$rfc);
  
  #todo Test we got the right wsize returned  
  #and the expected number of scores
  
  #Be careful we aren't just reimplementing the API method here
  #We could do this by subbing out more of the set_collection_config_by_Slice_ResultSets
  #method, and test each part separately here
  
  #my $max_bins_default = 700;
  #my $exp_wsize = ;
  my $wsize = $rfc->window_size;
  
  my $exp_scores = $slice->length / $rfc->window_size;
    
  if($exp_scores =~ /\./){ #we do not have an int
    $exp_scores = int($exp_scores) + 1; #int 'rounds' down and we want overlaps
  }
   
  ok(scalar(@{$rfc->scores}) == $exp_scores,
   "Found expected number($exp_scores) of scores for ResultFeature ".
   "Collection with window_size $wsize:\t".scalar(@{$rfc->scores}));   
  
}

  
#This actually needs  
my ($wsizes_valid, $wsizes);
  
SKIP:{
  skip "fetch_all_by_Slice_ResultSet valid ResultFeature Collection" if ! $rf_ref_valid;

  $wsizes       = $rf_a->window_sizes;
  $wsizes_valid = check_ref($wsizes, 'ARRAY');
  ok($wsizes_valid && (scalar(@$wsizes) == $num_window_sizes),
     "window_sizes returned ARRAYREF with expected number($num_window_sizes) of elements:".
     "\t$wsizes (".scalar(@$wsizes).')');
  
  #todo check they are all defined and ints
     
}
 
SKIP:{
  skip "Failed to get window_sizes" if ! $wsizes_valid;
    
  foreach my $wsize(@$wsizes){#8 wsize tests here
    #This is in ResultSet.t
    #my $col_file = $rset->get_dbfile_path_by_window_size($wsize);
    #ok(-f $col_file, "Collection (window_size=$wsize) file for $rset_name exists:\n".
    #  "\t$col_file");    
    
    #Now let's test the get_ResultFeatures_by_Slice wrapper for each window?
    ($rfcs) = values(%{$rf_a->fetch_all_by_Slice_ResultSets($slice, [$rset], undef, $wsize)});
    $rfc = $rfcs->[0];
  
    #don't have to test/skip the $rfc here as we have already done this outside the loop
   
    ok(($rfc->window_size == $wsize),
       "fetch_all_by_Slice_ResultSets - got expected window_size($wsize):\t".
       $rfc->window_size); 
       
    my $exp_scores = $slice->length / $wsize;
    #warn 'exp_scores '.$slice->length.' / '.$wsize.' = '.$exp_scores;
    
    if($exp_scores =~ /\./){ #we do not have an int
      $exp_scores = int($exp_scores) + 1; #int 'rounds' down and we want overlaps
    }
   
    #warn "exp_scores after round up $exp_scores";
    
    #This is currently returning 1 too many scores for the 65bp window size
    #which suggests there is an extra bin being erronously added somewhere
    #this appears not to be when the slice bounds match the bin location bounds???
   
 
    #warn "Slice:\t\t\t".$slice->start.' - '.$slice->end; #should be genomic
    #warn "ResultFeature:\t".$rfc->start.' - '.$rfc->end;  #should be local, no these are genomic!
    #todo change ResultFeatureAdaptor to return in local coords wrt query slice
   
   
    ok(scalar(@{$rfc->scores}) == $exp_scores,
       "Found expected number($exp_scores) of scores for ResultFeature ".
       "Collection with window_size $wsize:\t".scalar(@{$rfc->scores}));   
    
    #teststart is <= slice start
    
    #warn join(' ', @{$rfc->scores});
  }
  
  
  #Now do for fixed slice lengths and max bins to test the correct wsize used?
  
  #todo test the fetch_all_by_Slice_ResultSet wrapper method
  
}

done_testing();
