#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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
use autodie;

use Getopt::Long;
use Data::Printer;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw(throw);

use Bio::EnsEMBL::Funcgen::MirnaTargetFeature;
use Bio::EnsEMBL::DBEntry;

main();

sub main {

    my %options;
    GetOptions(\%options,
               'registry=s',
               'species=s',
               'input=s',
               'aliases=s',
               'release=s',
               'assembly=s'
    );

    my $miRBase_to_display_label = read_alias_file($options{'alias'});

    my $adaptors = get_adaptors(\%options);

    my $cache = {};
    $cache->{'feature_type'} =
        $adaptors->{'FeatureType'}->fetch_by_name('TarBase miRNA target');
    $cache->{'analysis'} =
        $adaptors->{'Analysis'}->fetch_by_logic_name('TarBase');

    open my $input_file, '<', $options{'input'};

    my $line_cnt = 0;
    my @features_to_store = ();
    my $batch_size = 10_000;
    while (readline $input_file) {
        my $line = $_;
        chomp $line;
        $line_cnt++;

        my $mirna_target_feature =
            parse_input_line($line, $line_cnt, \%options, $adaptors,
                             $miRBase_to_display_label, $cache);

        if($mirna_target_feature){
            push @features_to_store, $mirna_target_feature;
        }

        if (scalar @features_to_store == $batch_size){
            $adaptors->{'MirnaTargetFeature'}->store(@features_to_store);
            @features_to_store = ();
        }
    }
    $adaptors->{'MirnaTargetFeature'}->store(@features_to_store);

    close $input_file;
}

sub get_adaptors {
    my $options = shift;

    my $adaptors = {};

    Bio::EnsEMBL::Registry->load_all($options->{'registry'});

    my $funcgen_db =
        Bio::EnsEMBL::Registry->get_DBAdaptor($options->{'species'}, 'funcgen');
    my $core_db =
        Bio::EnsEMBL::Registry->get_DBAdaptor($options->{'species'}, 'core');

    my @adaptors_to_get = ('MirnaTargetFeature',
                           'Analysis',
                           'FeatureType'
    );

    for my $ad (@adaptors_to_get) {
        $adaptors->{$ad} = $funcgen_db->get_adaptor($ad);
    }

    $adaptors->{'gene'} = $core_db->get_adaptor('Gene');

    return $adaptors;
}

# MI0000027	cel-mir-56;
# MI0000028	cel-mir-57;
# MI0000029	cel-mir-58;cel-mir-58a;
sub read_alias_file {
    my ($file) = @_;

    my $miRBase_to_display_label = {};

    open my $fh, '<', $file;
    while (readline $fh) {
        chomp;
        my ($miRBase_accession, $labels_string) = split /\s+/, $_;
        my @labels                              = split /;/, $labels_string;

        # if there two or more, keep the last one
        my $display_label = pop @labels;

        $miRBase_to_display_label->{$miRBase_accession} = $display_label;
    }

    return $miRBase_to_display_label;
}

#0 MIMAT0000062|
#1 ENSG00000005810|
#2 HITS-CLIP|
#3 Experimental|
#4 77270375_77270394|
#5 BS2|
#6 spliced_no
#7 http://carolina.imis.athena-innovation.gr/diana_tools/web/index.php?r=tarbasev8%2Findex&miRNAs%5B%5D=MIMAT0000062&genes%5B%5D=ENSG00000005810
sub parse_input_line {
    my ($line, $line_cnt, $options, $adaptors, $miRBase_to_display_label,
        $cache) = @_;

    my @fields = split /\|/, $line;

    my $expected_fields = 8;
    my $actual_fields   = scalar @fields;
    if ($actual_fields != $expected_fields) {
        throw('Format error, expected ' . $expected_fields . ' fields, got '
            . $actual_fields . ' instead. Please check input line ' .
            $line_cnt . ' in ' . $options->{'input'});
    }

    my $mirna_target_feature =
        create_mirna_target_feature($adaptors,
                                    \@fields,
                                    $cache,
                                    $miRBase_to_display_label);

    return $mirna_target_feature;
}

sub create_mirna_target_feature {
    my ($adaptors, $fields, $cache, $miRBase_to_display_label) = @_;

    my $gene_stable_id = $fields->[1];
    if (!exists $cache->{'gene'}->{$gene_stable_id}) {
        $cache->{'gene'}->{$gene_stable_id} =
            $adaptors->{'gene'}->fetch_by_stable_id($gene_stable_id);
    }
    my $gene = $cache->{'gene'}->{$gene_stable_id};

    if (!$gene){
        warn('Gene ' . $gene_stable_id . ' not found!');
        return;
    }
    my $accession     = $fields->[0];
    my $display_label = $miRBase_to_display_label->{$accession};

    my ($start, $end) = split /_/, $fields->[4];

    my $mirna_target_feature =
        Bio::EnsEMBL::Funcgen::MirnaTargetFeature->new(
            -accession              => $accession,
            -method                 => $fields->[2],
            -evidence               => $fields->[3],
            -supporting_information => $fields->[5] . '; ' . $fields->[6],
            -feature_type           => $cache->{'feature_type'},
            -slice                  => $gene->slice(),
            -display_label          => $display_label,
            -start                  => $start,
            -end                    => $end,
            -strand                 => $gene->strand(),
            -analysis               => $cache->{'analysis'},
            -gene_stable_id         => $gene_stable_id,
        );

    return $mirna_target_feature;
}
