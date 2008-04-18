use strict;
use Bio::EnsEMBL::Registry;

my $reg = 'Bio::EnsEMBL::Registry';

$reg->load_registry_from_db(
							-host => 'ensembldb.ensembl.org',
							-user => 'anonymous',
						   );

my $efg_db = $reg->get_DBAdaptor('Human', 'funcgen');

#Investigate an Experiment
#Create a script which gets the 'ctcf_ren' Experiment.  Fetch all the ExperimentalChips which belong to this experiment.  List which feature and cell types were assayed using these chips.  Using the ExperimentalChip objects, retrieve the ArrayChip and the Array information corresponding to each ExperimentalChip.
#Hint: Try using the convinience methods in ExperimentalChip.

my $efg_db = $reg->get_DBAdaptor('Human', 'funcgen');

my $exp_adaptor = $efg_db->get_ExperimentAdaptor;

my $exp = $exp_adaptor->fetch_by_name('ctcf_ren');
my @exp_chips = @{$exp->get_ExperimentalChips};

print $exp->name.' '.$exp->primary_design_type.' experiment contains '.scalar(@exp_chips)." ExperimentalChips\n";

foreach my $exp_chip(@exp_chips){
  my $array_chip = $exp_chip->get_ArrayChip;
  my $array = $array_chip->get_Array;

  my $ctype = $exp_chip->cell_type;

  print 'ExperimentalChip '.$exp_chip->unique_id.' was used to assay '.$exp_chip->feature_type->name.' in '.$ctype->name.' '.$ctype->description."\n";
}
