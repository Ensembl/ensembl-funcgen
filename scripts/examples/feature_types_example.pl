use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

my $ftype_adaptor = $registry->get_adaptor('Human', 'funcgen', 'featuretype');
my $fset_adaptor = $registry->get_adaptor('Human', 'funcgen', 'featureset');

#Print all feature sets for Transcription Factors (note that this does not include CTCF, an insulator)
my @tfs = @{$ftype_adaptor->fetch_all_by_class('Transcription Factor')}; 
foreach my $ft (@tfs){
	print "Feature Type: ".$ft->name."\n";
	my @fsets = @{$fset_adaptor->fetch_all_by_FeatureType($ft)};
	print "\t".scalar(@fsets)." Feature Sets available:\n";
	foreach my $fset (@fsets){ 
		print "\t\t".$fset->name."\n"; 
	}
}
