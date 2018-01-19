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
    -user => 'anonymous'
);

# 2. Feature Types
# Get the list of all feature types in the Mouse eFG. For each one print its name, class and description
# Get the list of all Mouse feature types with class 'Transcription Factor'. How many transcription factors are there?

#Grab the eFG adaptor
my $fta = $registry->get_adaptor('Mouse', 'funcgen', 'featuretype');

my @mouse_feature_types = @{$fta->fetch_all()};

foreach my $feature_type (@mouse_feature_types){
	print $feature_type->name."\t".$feature_type->class."\t".$feature_type->description."\n";
}

my @tfs = @{$fta->fetch_all_by_class('Transcription Factor')};
print scalar(@tfs)." transcription factors in the mouse eFG\n";


__END__

>perl metadata_2.pl
H4ac    Histone Histone 4 Acetylation
H3ac    Histone Histone 3 Acetylation
H3K9ac  Histone Histone 3 Lysine 9 Acetylation
H3K4me3 Histone Histone 3 Lysine 4 Tri-Methylation
H3K4me2 Histone Histone 3 Lysine 4 Di-Methylation
H3K4me1 Histone Histone 3 Lysine 4 Mono-Methylation
H4K4me3 Histone Histone 4 Lysine 4 Tri-Methylation
Gene Associated Regulatory Feature      Gene like regulatory feature
Promoter Associated     Regulatory Feature      Promoter like regulatory feature
Non-Gene Associated     Regulatory Feature      Non-Gene like regulatory feature
Unclassified    Regulatory Feature      Unclassified regulatory feature
cisRED Search Region    Search Region   cisRED search region
cisRED  Regulatory Motif        cisRED group motif set
miRanda Target  RNA     miRanda microRNA target
VISTA Target    Search Region   VISTA target region
VISTA Enhancer  Enhancer        Enhancer identified by positive VISTA assay
VISTA Target - Negative Search Region   Enhancer negative region identified by VISTA assay
cisRED Motif    Regulatory Motif        cisRED motif
crtMmus200913   Regulatory Motif        cisRED group
crtMmus200731   Regulatory Motif        cisRED group
crtMmus200472   Regulatory Motif        cisRED group
crtMmus200738   Regulatory Motif        cisRED group
crtMmus200085   Regulatory Motif        cisRED group
crtMmus200761   Regulatory Motif        cisRED group
crtMmus200511   Regulatory Motif        cisRED group
crtMmus200727   Regulatory Motif        cisRED group
crtMmus200157   Regulatory Motif        cisRED group
crtMmus200777   Regulatory Motif        cisRED group
crtMmus200631   Regulatory Motif        cisRED group
crtMmus200042   Regulatory Motif        cisRED group
crtMmus201037   Regulatory Motif        cisRED group
crtMmus200736   Regulatory Motif        cisRED group
crtMmus200052   Regulatory Motif        cisRED group
crtMmus201023   Regulatory Motif        cisRED group
crtMmus200701   Regulatory Motif        cisRED group
crtMmus200394   Regulatory Motif        cisRED group
crtMmus200485   Regulatory Motif        cisRED group
crtMmus200792   Regulatory Motif        cisRED group
crtMmus200997   Regulatory Motif        cisRED group
crtMmus200450   Regulatory Motif        cisRED group
crtMmus200751   Regulatory Motif        cisRED group
crtMmus200773   Regulatory Motif        cisRED group
crtMmus200484   Regulatory Motif        cisRED group
crtMmus201009   Regulatory Motif        cisRED group
crtMmus200959   Regulatory Motif        cisRED group
crtMmus200449   Regulatory Motif        cisRED group
...
mmu-miR-872*    RNA     miRanda miRNA_target
mmu-miR-485*    RNA     miRanda miRNA_target
mmu-miR-294*    RNA     miRanda miRNA_target
mmu-miR-203*    RNA     miRanda miRNA_target
mmu-miR-196a*   RNA     miRanda miRNA_target
mmu-miR-879*    RNA     miRanda miRNA_target
mmu-miR-299*    RNA     miRanda miRNA_target
PolIII Transcription Associated Regulatory Feature      PolIII transcribed RNA gene associated regulatory feature
cMyb    Transcription Factor    cMyb Transcription Factor Binding
Max     Transcription Factor    Max Transcription Factor Binding
NELFe   Transcription Factor    NELFe Transcription Factor Binding
Rad21   Transcription Factor    Rad21 Transcription Factor Binding
USF2    Transcription Factor    USF2 Transcription Factor Binding
GFP     DNA     GFP Control
Rbbp5   Transcription Factor    Rbbp5 Transcription Factor Binding
Wdr5    Transcription Factor    Wdr5 Transcription Factor Binding
negative        DNA     negative library

23 transcription factors in the mouse eFG DB
