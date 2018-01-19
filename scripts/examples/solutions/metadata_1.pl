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

# 1. Cell Types
# Get the list of all cell types in the Human eFG. For each one print its name, gender and description
# Hint: not all cells have a determined gender

#Grab the eFG adaptor
my $cta = $registry->get_adaptor('Human', 'funcgen', 'celltype');

my @human_cell_types = @{$cta->fetch_all()};

foreach my $cell (@human_cell_types){
	my $gender = $cell->gender || "Undetermined gender"; 
	print $cell->name."\t".$gender."\t".$cell->description."\n";
}

__END__

>perl metadata_1.pl

HeLa-S3 female  Human Epithelial Carcinoma Cells
GM06990 female  Human B-Lymphocyte cell line
U2OS    female  Human Bone Osteosarcoma Epithelial Cells
CD4     Undetermined gender     Human CD4 T-Cells
IMR90   female  Human Fetal Lung Fibroblast
HL-60   female  Human Promyelotic Leukemia Cells
HepG2   male    Human hepatocellular liver carcinoma cell line
Lymphoblastoid  Undetermined gender     Lymphoblastoid cells
CD133   Undetermined gender     Human CD133+ Hematopoietic Stem Cell
CD36    Undetermined gender     Human CD36+ Erythrocyte Precursor Cell
K562    female  Human myelogenous leukaemia cell line
GM12878 female  Human B-Lymphocyte cell line
HUVEC   Undetermined gender     Human umbilical vein endothelial cell line
NHEK    female  Normal human epidermal keratinocyte cell line
H1ESC   Undetermined gender     Human Embryonic Stem Cells
MultiCell       Undetermined gender     Multiple CellTypes used in core RegulatoryFeature set
K562b   female  Human myelogenous leukaemia line (alternative)
NH-A    Undetermined gender     Normal Human Astrocytes
HSMM    Undetermined gender     Human Skeletal Muscle Myoblasts
HMEC    Undetermined gender     Human Mammary Epithelial Cells
