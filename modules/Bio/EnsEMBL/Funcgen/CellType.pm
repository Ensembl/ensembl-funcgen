#
# Ensembl module for Bio::EnsEMBL::Funcgen::CellType
#


=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

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


=head1 NAME

Bio::EnsEMBL::Funcgen::CellType - A module to represent a CellType.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::CellType;

#Fetch from adaptor
my $ctype = $cell_type_adaptor->fetch_by_name($ctype_name);

#Create from new
my $ctype = Bio::EnsEMBL::Funcgen::CellType->new
  (
   -name          => 'H1-ESC',
   -display_label => 'H1-ESC',
   -description   => 'Human Embryonic Stem Cell',
   -efo_id        => 'efo:EFO_0003042',
   -tissue        => 'embryonic stem cell',
  );

print $ctype->name.' is a '.$ctype->description."\n";

#H1-ESC is a Human Embryonic Stem Cell


=head1 DESCRIPTION

This is a simple class to represent information about a cell type.  This may represent
harvested cells, a cell line or a more generic tissue type.

=cut

package Bio::EnsEMBL::Funcgen::CellType;

use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw( deprecate );
use base qw(Bio::EnsEMBL::Funcgen::Epigenome);

deprecate(
    "Bio::EnsEMBL::Funcgen::CellType has been deprecated and will be removed in Ensembl release 89."
        . " Please use Bio::EnsEMBL::Funcgen::Epigenome instead."
);

1;

