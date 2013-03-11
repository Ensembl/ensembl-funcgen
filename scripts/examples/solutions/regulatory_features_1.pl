#!/usr/bin/perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
);

# Regulatory Features 
# Using the human DB, fetch the all the cell type specific regulatory features with stable ID 'ENSR00000623613'.
# Print out the stable ID, bound_start/end and start/end values, name of the cell and feature type for each.
# HINT: To get all the cell type specific RegulatoryFeatures use the fetch_all_by_stable_ID method.

#Grab the regfeat adaptor
my $regfeat_adaptor = $registry->get_adaptor('Human', 'funcgen', 'regulatoryfeature');


#Fetch ENSR00000623613 reg feats for all cell types.
my @rfs = @{$regfeat_adaptor->fetch_all_by_stable_ID('ENSR00000623613')};


# Print out the details
foreach my $rf(@rfs){
  print $rf->stable_id.": \n";
  print "\tRelative Position: ".$rf->bound_start."..".$rf->start."-".$rf->end."..".$rf->bound_end."\n";
  print "\tCell: ".$rf->cell_type->name."\n";
  print "\tFeature Type: ".$rf->feature_type->name."\n";
}

__END__


ENSR00000623613: 
        Relative Position: 27178100..27178438-27185061..27186300
        Cell: K562
        Feature Type: Unclassified
ENSR00000623613: 
        Relative Position: 27123400..27178438-27185061..27297450
        Cell: H1ESC
        Feature Type: Unclassified
ENSR00000623613: 
        Relative Position: 27177575..27178438-27185061..27187569
        Cell: HUVEC
        Feature Type: Gene Associated
ENSR00000623613: 
        Relative Position: 27177907..27178438-27185061..27185882
        Cell: HMEC
        Feature Type: Promoter Associated
ENSR00000623613: 
        Relative Position: 27126300..27178438-27185061..27270800
        Cell: NH-A
        Feature Type: Unclassified
ENSR00000623613: 
        Relative Position: 27178438..27178438-27185061..27185061
        Cell: HeLa-S3
        Feature Type: Unclassified
ENSR00000623613: 
        Relative Position: 27178438..27178438-27185061..27187350
        Cell: HepG2
        Feature Type: Unclassified
ENSR00000623613: 
        Relative Position: 27178050..27178438-27185061..27187700
        Cell: IMR90
        Feature Type: Unclassified
ENSR00000623613: 
        Relative Position: 27178438..27178438-27185061..27185061
        Cell: HSMM
        Feature Type: Promoter Associated
ENSR00000623613: 
        Relative Position: 27178438..27178438-27185061..27185061
        Cell: MultiCell
        Feature Type: Unclassified
ENSR00000623613: 
        Relative Position: 27177500..27178438-27185061..27187050
        Cell: NHEK
        Feature Type: Unclassified
ENSR00000623613: 
        Relative Position: 27177650..27178438-27185061..27185061
        Cell: GM12878
        Feature Type: Unclassified
