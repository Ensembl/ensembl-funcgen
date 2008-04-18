use strict;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db
  (
   -host => 'ensembldb.ensembl.org',
   -user => 'anonymous',
  );


#Grab the adaptors
my $efg_db           = $registry->get_DBAdaptor('Human', 'funcgen');
my $exp_adaptor      = $efg_db->get_ExperimentAdaptor;

#Grab a the CTCF experiment
my $exp              = $exp_adaptor->fetch_by_name('ctcf_ren');

#Print some info
my $num_chips = scalar(@{$exp->get_ExperimentalChips});
print $exp->name.' '.$exp->primary_design_type." experiment contains $num_chips ExperimentalChips\n";




#Grab the ResultSetAdaptor
my $resultset_adaptor = $efg_db->get_ResultSetAdaptor;
my $slice_adaptor     = $efg_db->get_SliceAdaptor;

#Grab a slice for chr X
my $slice = $slice_adaptor->fetch_by_region('chromosome', 'X');

#Grab the ctcf biological replicate 1 result set
#We are assigning the array to one element as we know there is 
#only one entry in the DB for this name
my ($result_set) = @{$resultset_adaptor->fetch_all_by_name('ctcf_ren_BR1')};

my @result_features = @{$result_set->get_ResultFeatures_by_Slice($slice)};
print "Chromosome X has ".scalar(@result_features)." results\n";

my $count = 0;

foreach my $rf(@result_features){
  print "Locus:\t".$rf->start.'-'.$rf->end."\tScore:".$rf->score."\n";

  $count ++;
  last if $count == 10;
}

