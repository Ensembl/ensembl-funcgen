#
# Ensembl module for Bio::EnsEMBL::Funcgen::BindingMatrix::Constants
#

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Funcgen::DBSQL::BindingMatrixConstants

=head1 DESCRIPTION

The BindingMatrixConstants is a module for declaring the unit constants of the
BindingMatrix class


=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::BindingMatrix
Bio::EnsEMBL::Funcgen::BindingMatrix::Converter

=cut

package Bio::EnsEMBL::Funcgen::BindingMatrix::Constants;

use base qw( Exporter );
use vars qw( @EXPORT_OK );
use strict;
use warnings;

our @EXPORT_OK = qw(
    FREQUENCIES
    PROBABILITIES
    WEIGHTS
    BITS
    VALID_UNITS
    VALID_NUCLEOTIDES
);

our %EXPORT_TAGS = ( all => \@EXPORT_OK );

use constant {
    FREQUENCIES       => 'Frequencies',
    PROBABILITIES     => 'Probabilities',
    WEIGHTS           => 'Weights',
    BITS              => 'Bits',
    VALID_UNITS       => [ 'Frequencies', 'Probabilities', 'Weights', 'Bits' ],
    VALID_NUCLEOTIDES => [ 'A', 'C', 'G', 'T' ]
};

1;
