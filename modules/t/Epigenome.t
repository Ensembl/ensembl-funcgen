# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
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

# use Bio::EnsEMBL::Test::MultiTestDB;

# Module compiles
BEGIN { use_ok('Bio::EnsEMBL::Funcgen::Epigenome'); }

# Test constructor
my $epigenome = Bio::EnsEMBL::Funcgen::Epigenome->new(
    -name               => 'H1-ESC',
    -display_label      => 'H1-ESC',
    -description        => 'Human Embryonic Stem Cell',
    -gender             => 'female',
    -ontology_accession => 'efo:EFO_0003042',
    -tissue             => 'embryonic stem cell',
);

isa_ok( $epigenome, 'Bio::EnsEMBL::Funcgen::Epigenome', 'Epigenome' );

# Test name definition
throws_ok { Bio::EnsEMBL::Funcgen::Epigenome->new }
qr/Must supply an Epigenome name/, 'Check that name is supplied';

# Test gender definition
throws_ok {
    Bio::EnsEMBL::Funcgen::Epigenome->new(
        -name               => 'H1-ESC',
        -display_label      => 'H1-ESC',
        -description        => 'Human Embryonic Stem Cell',
        -gender             => 'invalid',
        -ontology_accession => 'efo:EFO_0003042',
        -tissue             => 'embryonic stem cell',
    );
}
qr/Gender .+ not valid, must be one of/, 'Check that the gender is valid';

# Test getter subroutines
is( $epigenome->name,          'H1-ESC', 'Retrieve epigenome name' );
is( $epigenome->display_label, 'H1-ESC', 'Retrieve epigenome display name' );
is( $epigenome->description,
    'Human Embryonic Stem Cell',
    'Retrieve epigenome description'
);
is( $epigenome->gender, 'female', 'Retrieve epigenome gender' );
is( $epigenome->ontology_accession,
    'efo:EFO_0003042', 'Retrieve epigenome ontology_accession' );
is( $epigenome->tissue, 'embryonic stem cell', 'Retrieve epigenome tissue' );

done_testing();
