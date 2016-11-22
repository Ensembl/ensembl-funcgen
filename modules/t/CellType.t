# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

use strict;
use warnings;
use Test::More;
use Test::Exception;
use Bio::EnsEMBL::Test::TestUtils;
use Bio::EnsEMBL::Funcgen::CellType;

# use Bio::EnsEMBL::Test::MultiTestDB;

# Module compiles
BEGIN { use_ok('Bio::EnsEMBL::Funcgen::CellType'); }

# Test constructor
my $ct = Bio::EnsEMBL::Funcgen::CellType->new(
    -name               => 'H1-ESC',
    -display_label      => 'H1-ESC',
    -description        => 'Human Embryonic Stem Cell',
    -gender             => 'female',
    -ontology_accession => 'efo:EFO_0003042',
    -tissue             => 'embryonic stem cell',
);

isa_ok( $ct, 'Bio::EnsEMBL::Funcgen::CellType', 'CellType' );

# Test name definition
throws_ok { Bio::EnsEMBL::Funcgen::CellType->new }
qr/Must supply an Epigenome name/, 'Check that name is supplied';

# Test gender definition
#throws_ok {
#    Bio::EnsEMBL::Funcgen::CellType->new(
#        -name               => 'H1-ESC',
#        -display_label      => 'H1-ESC',
#        -description        => 'Human Embryonic Stem Cell',
#        -gender             => 'invalid',
#        -ontology_accession => 'efo:EFO_0003042',
#        -tissue             => 'embryonic stem cell',
#    );
#}
#qr/Gender not valid, must be one of/, 'Check that the gender is valid';

# Test getter subroutines
is( $ct->name,          'H1-ESC', 'Retrieve cell type name' );
is( $ct->display_label, 'H1-ESC', 'Retrieve cell type display name' );
is( $ct->description,
    'Human Embryonic Stem Cell',
    'Retrieve cell type description'
);
is( $ct->gender, 'female', 'Retrieve cell type gender' );
is( $ct->ontology_accession, 'efo:EFO_0003042',
    'Retrieve cell type ontology_accession' );
is( $ct->tissue, 'embryonic stem cell', 'Retrieve cell type tissue' );

done_testing();
