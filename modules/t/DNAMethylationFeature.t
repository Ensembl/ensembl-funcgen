# Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

#!usr/bin/env perl
use strict;
use warnings;
use Test::More;
use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::Funcgen::DNAMethylationFeature;

#obtain Adaptors for dnabb and funcgen databases
my $multi  = Bio::EnsEMBL::Test::MultiTestDB->new();
my $efgdba = $multi->get_DBAdaptor("funcgen");

#define data root
my $data_root = '';#NEEDS DEFINING

throw('$data_root needs defining');

#obtain Adaptors
my $rsa           = $efgdba->get_adaptor("resultset");
my $dnaa          = $efgdba->get_adaptor("DNAMethylationFeature");
my $slice_adaptor = $efgdba->get_adaptor("slice");

#set data root
$rsa->dbfile_data_root($data_root);

#Get all ResultSets from DB that refer to DNAMethylationFeatures

my @result_sets = @{
    $rsa->fetch_all_by_feature_class( 'dna_methylation',
        { status => 'DISPLAYABLE' } )
};

plan tests => ( scalar @result_sets * 8 ) + 8;

my $slice = $slice_adaptor->fetch_by_region( 'chromosome', 1, 1, 3011550 );

foreach my $resultset (@result_sets) {

    my $constraints = {
        min_read_depth  => 15,
        context         => 'CG',      #'CT',
        min_methylation => '80',      #-1    # Float percentage
        max_methylation => '100',    #110 # Float percentage
    };

    #Need to do with and without constraints
    #as without uses a different constructor wrapper method

    my $dna_meth_features =
      $dnaa->fetch_all_by_Slice_ResultSet( $slice, $resultset, $constraints );

    ok( @{$dna_meth_features}, 'fetch_all_by_Slice_ResultSet' );

    #Test one random DNAMethylationFeature object
    my $index = int( rand( scalar @{$dna_meth_features} ) );

    my $df = $dna_meth_features->[$index];

    ok( $df->isa('Bio::EnsEMBL::Funcgen::DNAMethylationFeature'),
        "class - DNAMethylationFeature" );

    #This is not quite true as there maybe many other DNA meth analyses

    ok(
        (
                 $df->analysis->logic_name eq 'RRBS_merged_filtered_10'
              or $df->analysis->logic_name eq 'Whole_Genome_methylC_seq'
              or $df->analysis->logic_name eq 'WGBS_merged'
        ),
        "analysis logic_name: " . $df->analysis->logic_name
    );

    # ok( $df->context =~ m/CG|CHG|CHH/, "sequence context - " . $df->context );

    is(
        $df->methylated_reads,
        $df->total_reads - $df->unmethylated_reads,
        'read numbers consistancy'
    );

    cmp_ok( $df->total_reads, '>=', $constraints->{min_read_depth},
        'min_read_depth' );
    cmp_ok( $df->context, 'eq', $constraints->{context}, 'context' );
    cmp_ok( $df->percent_methylation, '>=', $constraints->{min_methylation},
        'min_methylation' );
    cmp_ok( $df->percent_methylation, '<=', $constraints->{max_methylation},
        'max_methylation' );

}

#Create and test a new DNAMethylationFeature object
my $dnamethylationfeature = Bio::EnsEMBL::Funcgen::DNAMethylationFeature->new(
    -SLICE               => $slice,
    -STRAND              => 1,
    -START               => 1,
    -END                 => 2,
    -METHYLATED_READS    => 33,
    -TOTAL_READS         => 37,
    -PERCENT_METHYLATION => 89,
    -DISPLAY_LABEL       => "display label test",
    -CONTEXT             => "CG",
    -ADAPTOR             => $dnaa,
    -SET                 => $result_sets[0],
);

ok( $dnamethylationfeature, 'new DNAMethylationFeature' );

ok( $dnamethylationfeature->isa('Bio::EnsEMBL::Funcgen::DNAMethylationFeature'),
    "class - DNAMethylationFeature" );

is( $dnamethylationfeature->methylated_reads, 33, 'methylated_reads' );

is( $dnamethylationfeature->unmethylated_reads, 4, 'unmethylated_reads' );

is( $dnamethylationfeature->total_reads, 37, 'total_reads' );

is( $dnamethylationfeature->percent_methylation, 89, 'percent_methylation' );

is( $dnamethylationfeature->context, 'CG', 'context' );

is(
    $dnamethylationfeature->display_label,
    'display label test',
    'display_label'
);
