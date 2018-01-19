#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
);

# Regulatory Feature vs ENCODE Segmentation classification


#Grab the adaptors
my $regfeat_adaptor = $registry->get_adaptor('Human', 'funcgen', 'regulatoryfeature');
my $segfeat_adaptor = $registry->get_adaptor('Human', 'funcgen', 'segmentationfeature');
my %ctype_segfeats  = ();



#Fetch ENSR00000623613 reg feats for all cell types.
my @rfs = @{$regfeat_adaptor->fetch_all_by_stable_ID('ENSR00001348194')};



print "Regulatory Build vs ENCODE segmentation classifications for:\t".
  $rfs[0]->stable_id.' ('.$rfs[0]->length."bp)\n";
my $rf_slice = $rfs[0]->feature_Slice;

foreach my $seg_feat( @{$segfeat_adaptor->fetch_all_by_Slice($rf_slice)} ){
  $ctype_segfeats{ $seg_feat->cell_type->name } ||= [];

  push @{ $ctype_segfeats{ $seg_feat->cell_type->name } }, $seg_feat;
}



# Print out the details
foreach my $rf(@rfs){
  my $ctype = $rf->cell_type->name;

  print "\nCellType\t".$ctype."\n";
  print "RegulatoryFeature Classification:\t".$rf->feature_type->name."\n";


  if(exists $ctype_segfeats{$ctype}){

    
    print "SegmentationFeature Classification:\t".
      join(', ', ( map { $_->feature_type->name } 
                   @{ $ctype_segfeats{ $ctype } } 
                 ))."\n";
  }
  else{
    print "SegmentationFeature Classication:\tNo ENCODE segmentation available for ".$ctype."\n";
  }
}



__END__


perl regulatory_features_4.pl
Regulatory Build vs ENCODE segmentation classifications for:    ENSR00001348194 (1530bp)

CellType        NHEK
RegulatoryFeature Classification:       Promoter Associated
SegmentationFeature Classication:       No ENCODE segmentation available for NHEK

CellType        H1ESC
RegulatoryFeature Classification:       Promoter Associated
SegmentationFeature Classication:       Predicted Repressed/Low Activity, Predicted Enhancer, Predicted Repressed/Low Activity

CellType        HUVEC
RegulatoryFeature Classification:       Gene Associated
SegmentationFeature Classication:       Predicted Promoter with TSS

CellType        HMEC
RegulatoryFeature Classification:       Promoter Associated
SegmentationFeature Classication:       No ENCODE segmentation available for HMEC

CellType        CD4
RegulatoryFeature Classification:       Promoter Associated
SegmentationFeature Classication:       No ENCODE segmentation available for CD4

CellType        NH-A
RegulatoryFeature Classification:       Gene Associated
SegmentationFeature Classication:       No ENCODE segmentation available for NH-A

CellType        HepG2
RegulatoryFeature Classification:       Unclassified
SegmentationFeature Classication:       Predicted Repressed/Low Activity

CellType        IMR90
RegulatoryFeature Classification:       Unclassified
SegmentationFeature Classication:       No ENCODE segmentation available for IMR90

CellType        HSMM
RegulatoryFeature Classification:       Promoter Associated
SegmentationFeature Classication:       No ENCODE segmentation available for HSMM

CellType        MultiCell
RegulatoryFeature Classification:       Unclassified
SegmentationFeature Classication:       No ENCODE segmentation available for MultiCell

CellType        K562
RegulatoryFeature Classification:       Unclassified
SegmentationFeature Classication:       Predicted Promoter with TSS

CellType        GM12878
RegulatoryFeature Classification:       Promoter Associated
SegmentationFeature Classication:       Predicted Promoter with TSS


#http://www.ensembl.org/Homo_sapiens/Share/c4f593fa0b3bb188eb088714440f380b89648374



#Ensembl Regulatory Build lacks specificity of classification with many unclassified, also has potential resolution issues. Also does no explicitly annotate 'dead' regions/
#Conversely...
#ENCODE segmentation lacks restricted to 6 cell types, as well as potential fragmention issues. No common feature location definition across cell lines. No MultiCell summary set. ~Genomewide coverage,


