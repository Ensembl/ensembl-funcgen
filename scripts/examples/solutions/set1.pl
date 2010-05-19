use strict;
use Bio::EnsEMBL::Registry;

my $reg = 'Bio::EnsEMBL::Registry';

$reg->load_registry_from_db(
							-host => 'ensembldb.ensembl.org',
							-user => 'anonymous',
						   );

my $efg_db = $reg->get_DBAdaptor('Human', 'funcgen');

#Investigate what DataSets are present in the data base.
#Create a script which lists all the DataSet names, feature types, cell types and their supporting set type



my $dataset_adaptor = $efg_db->get_DataSetAdaptor;

my @data_sets = @{$dataset_adaptor->fetch_all};

foreach my $dset(@data_sets){
  my $cell_info;

  my $ctype = $dset->cell_type;
  
  if($ctype){
    $cell_info = ' in '.$ctype->name.' '.$ctype->description;
  }

  print $dset->name.' contains '.$dset->feature_type->name.' features'.$cell_info."\n";
  #' from supporting set type '.ucfirst($dset->supporting_set_type)."Set\n";

} 
