#!/usr/bin/perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
);

# Regulatory Feature vs ENCODE Segmentation classification


#Lack of specific classification is an issue for the RegulatoryBuild, as well as potential resolution issues.
#Conversely...
#Cell type coverage is an issue for the ENCODE segmentation, as well as potential fragmention issues.

#Grab the adaptors
my $regfeat_adaptor    = $registry->get_adaptor('Human', 'funcgen', 'regulatoryfeature');
my $featureset_adaptor = $registry->get_adaptor('Human', 'funcgen', 'featureset');
my %segfeat_fsets = ();


# Cache SegmentatioFeature FeatureSets
foreach my $fset( @{$featureset_adaptor->fetch_all_by_feature_class('segmentation')} ){
  $segfeat_fsets{$fset->cell_type->name} = $fset;
}


#Fetch ENSR00000623613 reg feats for all cell types.
my @rfs = @{$regfeat_adaptor->fetch_all_by_stable_ID('ENSR00000623613')};



print "Regulatory Build vs ENCODE segmentation classifications for:\t".
  $rfs[0]->stable_id.' ('.$rfs[0]->length."bp)\n";
my $rf_slice = $rfs[0]->feature_Slice;

# Print out the details
foreach my $rf(@rfs){
  print "\nCellType\t".$rf->cell_type->name."\n";
  print "RegulatoryFeature Classification:\t".$rf->feature_type->name."\n";


  if(exists $segfeat_fsets{$rf->cell_type->name}){

    my @seg_feats = @{$segfeat_fsets{$rf->cell_type->name}->
                        get_Features_by_Slice($rf_slice)};
    
    print "SegmentationFeature Classication:\t".
      join(', ', ( map { $_->feature_type->name } @seg_feats))."\n";
  }
  else{
    print "SegmentationFeature Classication:\tNo ENCODE segmentation available for ".$rf->cell_type->name."\n";
  }


}

__END__
