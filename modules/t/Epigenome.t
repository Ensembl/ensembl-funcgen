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
use Bio::EnsEMBL::Test::MultiTestDB;

# Setup connection to test database
my $multi = Bio::EnsEMBL::Test::MultiTestDB->new();
my $funcgen_db = $multi->get_DBAdaptor('funcgen');

# Module compiles
BEGIN { use_ok('Bio::EnsEMBL::Funcgen::Epigenome'); }

# Test constructor
my $new_epigenome = Bio::EnsEMBL::Funcgen::Epigenome->new(
    -name               => 'H1-ESC',
    -display_label      => 'H1-ESC',
    -description        => 'Human Embryonic Stem Cell',
    -gender             => 'female',
    -ontology_accession => 'efo:EFO_0003042',
    -tissue             => 'embryonic stem cell',
);

isa_ok( $new_epigenome, 'Bio::EnsEMBL::Funcgen::Epigenome', 'Epigenome' );

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
my $epigenome_adaptor = $funcgen_db->get_adaptor('Epigenome');
my $epigenome = $epigenome_adaptor->fetch_by_name('NHDF-AD');

is( $epigenome->name,          'NHDF-AD', 'Test name()' );
is( $epigenome->production_name,          'NHDF_AD', 'Test production_ame()' );
is( $epigenome->display_label, 'NHDF-Ad', 'Test display_label()' );
is( $epigenome->description,
    'adult dermal fibroblasts',
    'Test description()'
);
is( $epigenome->gender, 'female', 'Test gender()' );
is($epigenome->efo_accession, undef, 'Test efo_accession()');
is($epigenome->epirr_accession, undef, 'Test epirr_accession()');

my $expected_summary = {
    name          => 'NHDF-AD',
    gender        => 'female',
    description   => 'adult dermal fibroblasts',
    display_label => 'NHDF-Ad',
    efo_accession => undef
};

is_deeply( $epigenome->summary_as_hash, $expected_summary,
    'Test summary_as_hash()' );

done_testing();
