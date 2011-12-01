#!/usr/bin/perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

# Annotated Features
# Annotated Features result from analysis of raw data, corresponding to regions in the genome enriched for specific events (like TF binding or Histone Marks)

#Grab the eFG adaptor
my $sa = $registry->get_adaptor('Human', 'core', 'slice');
my $fsa = $registry->get_adaptor('Human', 'funcgen', 'featureset');

# You can get it by the name
#Carefull: some parts of the Y chomosome are the same as the Y...
my $slice = $sa->fetch_by_region('chromosome','Y',5000000,40000000);

# Compare the number of annotated features in the region Y:5000000-40000000 between the Human feature sets 'K562_DNase1_ENCODE_Duke_SWEmbl_R0025_D150' and 'HepG2_DNase1_ENCODE_Duke_SWEmbl_R0025_D150'.Investigate the differences (hint: check cell types).

my $K562_DNAse_fset = $fsa->fetch_by_name('K562_DNase1_ENCODE_Duke_SWEmbl_R0025_D150');
my @K562_DNAse_afs = @{$K562_DNAse_fset->get_Features_by_Slice($slice)};
print "K562 contains ".scalar(@K562_DNAse_afs)." features\n";
print $K562_DNAse_fset->cell_type->name." has gender ".$K562_DNAse_fset->cell_type->gender."\n";
my $HepG2_DNAse_fset = $fsa->fetch_by_name('HepG2_DNase1_ENCODE_Duke_SWEmbl_R0025_D150');
my @HepG2_DNAse_afs = @{$HepG2_DNAse_fset->get_Features_by_Slice($slice)};
print "HepG2 contains ".scalar(@HepG2_DNAse_afs)." features\n";
print $HepG2_DNAse_fset->cell_type->name." has gender ".$HepG2_DNAse_fset->cell_type->gender."\n";

# Print the position, score and summit of annotated features for feature set 'K562_CTCF_ENCODE_Broad_SWEmbl_R015_D150' within the region 1:100000-150000
print "Features for K562_CTCF_ENCODE_Broad_SWEmbl_R015_D150\n";
my $K562_CTCF_fset = $fsa->fetch_by_name('K562_CTCF_ENCODE_Broad_SWEmbl_R015_D150');
$slice = $sa->fetch_by_region('chromosome','1',100000,150000);
my @K562_CTCF_afs = @{$K562_CTCF_fset->get_Features_by_Slice($slice)};
foreach my $af (@K562_CTCF_afs){
	print "Start: ".$af->seq_region_start." End: ".$af->seq_region_end." Score: ".$af->score." Summit: ".$af->summit."\n";
}
print "\n";



# (optional) Print the properties of the analysis used in the previous feature sets
print "K562 DNAse\n";
print_analysis($K562_DNAse_fset->analysis);
print "\n";

print "HepG2 DNAse\n";
print_analysis($HepG2_DNAse_fset->analysis);
print "\n";

print "K562 CTCF\n";
print_analysis($K562_CTCF_fset->analysis);
print "\n";

sub print_analysis {
	my $an = shift;
	print "Logic Name: ".$an->logic_name."\n";
	print "Display Label: ".$an->display_label."\n";
	print "Program: ".$an->program."\t".$an->parameters."\n";
}
