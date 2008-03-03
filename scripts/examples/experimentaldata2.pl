use strict;
use Bio::EnsEMBL::Registry;

my $reg = 'Bio::EnsEMBL::Registry';

$reg->load_registry_from_db(
							-host => 'ensembldb.ensembl.org',
							-user => 'anonymous',
						   );

my $efg_db = $reg->get_DBAdaptor('Human', 'funcgen');


#Investigate the ResultSets for this Experiment
#Create a script which retrieves all the ResultSets for the 'ctcf_ren' Experiment, list their names and their analysis logic_names.  Using the normalised ResultSet named 'ctcf_ren_BR1_TR1', retrieve all the ResultFeatures for the slice 'chromosome:NCBI36:X:27126000:27169000'.  Print out how many ResultFeatures were returned and list the start, end and score values for the first 5.

my $efg_db = $reg->get_DBAdaptor('Human', 'funcgen');
my $exp_adaptor = $efg_db->get_ExperimentAdaptor;
my $resultset_adaptor = $efg_db->get_ResultSetAdaptor;

my $exp = $exp_adaptor->fetch_by_name('ctcf_ren');
my @rsets = @{$resultset_adaptor->fetch_all_by_Experiment($exp)};


foreach my $rset(@rsets){
  print 'ResultSet '.$rset->name.' was analysed using '.$rset->analysis->logic_name."\n";
}

my ($rset) = @{$resultset_adaptor->fetch_all_by_name('ctcf_ren_BR1_TR1')};
#There can be more than one as we might have different analyses for the same set of data.
#In this case, there is only one.

my $slice = $efg_db->get_SliceAdaptor->fetch_by_name('chromosome:NCBI36:X:27126000:27169000');
my @result_features = @{$rset->get_ResultFeatures_by_Slice($slice)};
print "Found ".scalar(@result_features).' ResultFeatures for slice '.$slice->name."\n";

my $cnt =0;

foreach my $rfeature(@result_features){

  print 'Start: '.$rfeature->start.' End: '.$rfeature->end.' Score: '.$rfeature->score."\n";

  $cnt++;
  last if $cnt ==5;
}

