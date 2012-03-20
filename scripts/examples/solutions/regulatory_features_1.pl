#!/usr/bin/perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

# Regulatory Features 
# Create a script which gets all the RegulatoryFeatures within 2KB of the 'HOXA4' gene.
# Hint: Use the Gene->feature_Slice method, then use the expand method.
# Print out their stable IDs, bound_start/end and start/end values, name of the cell and feature types.

#Grab the eFG adaptor
my $gene_adaptor = $registry->get_adaptor('Human', 'core', 'gene');
my $regfeat_adaptor = $registry->get_adaptor('Human', 'funcgen','regulatoryfeature');

my @genes = @{$gene_adaptor->fetch_all_by_external_name('HOXA4')};
print scalar(@genes)." human gene(s) named Hox4\n";

my $hox4 = $genes[0];
my $slice = $hox4->feature_Slice;
print "Slice corresponding to Hox4: ".$slice->seq_region_name."\t".$slice->start."\t".$slice->end."\t".$slice->strand."\n";
$slice = $slice->expand(2000,2000);
print "Slice 2Kb around Hox4: ".$slice->seq_region_name."\t".$slice->start."\t".$slice->end."\t".$slice->strand."\n";

my @reg_feats = @{$regfeat_adaptor->fetch_all_by_Slice($slice)};
print scalar(@reg_feats)." Regulatory Features 2kb around Hox4\n";

foreach my $rf (@reg_feats){
	print $rf->stable_id.": \n";
	print "\tRelative Position: ".$rf->bound_start."..".$rf->start."-".$rf->end."..".$rf->bound_end."\n";
	print "\tCell: ".$rf->cell_type->name."\n";
	print "\tFeature Type: ".$rf->feature_type->name."\n";
}


__END__


>perl regulatory_features_1.pl

1 human gene(s) named Hox4
Slice corresponding to Hox4: 7  27168126        27170418        -1
Slice 2Kb around Hox4: 7        27166126        27172418        -1

4 Regulatory Features 2kb around Hox4
ENSR00001558218:
        Relative Position: 5695..5695-5950..5950
        Cell: MultiCell
        Feature Type: Unclassified
ENSR00001049354:
        Relative Position: 4335..4335-4845..4845
        Cell: MultiCell
        Feature Type: Unclassified
ENSR00000623608:
        Relative Position: 1188..1188-2849..2849
        Cell: MultiCell
        Feature Type: Unclassified
ENSR00001049355:
        Relative Position: 296..296-560..560
        Cell: MultiCell
        Feature Type: Unclassified
