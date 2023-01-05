#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2023] EMBL-European Bioinformatics Institute

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
use Bio::EnsEMBL::Utils::SqlHelper;

use Bio::EnsEMBL::Funcgen::MirnaTargetFeature;
use Bio::EnsEMBL::DBEntry;

main();

sub main {

    my %options;
    GetOptions(\%options,
               'registry=s',
               'species=s',
               'input=s',
               'assembly=s'
    );

    my $adaptors = get_adaptors(\%options);

    my $cache = {};

    $cache->{feature_set} =
        $adaptors->{'FeatureSet'}->fetch_by_name('VISTA enhancer set');

    my $helper = Bio::EnsEMBL::Utils::SqlHelper->new(
        -DB_CONNECTION => $adaptors->{'ExternalFeature'}->db->dbc
    );

    my $sql_command = 'DELETE FROM external_feature WHERE feature_set_id='
        . $cache->{feature_set}->dbID();
    $helper->execute_simple(-SQL => $sql_command);

    $cache->{'positive_feature_type'} =
        $adaptors->{'FeatureType'}->fetch_by_name('VISTA Enhancer');
    $cache->{'negative_feature_type'} =
        $adaptors->{'FeatureType'}->fetch_by_name('VISTA Target - Negative');

    open my $input_file, '<', $options{'input'};

    my @features_to_store = ();
    my $batch_size        = 10_000;

    while (readline $input_file) {
        my $line = $_;
        chomp $line;

        my $vista_feature =
            parse_input_line($line, $adaptors, $cache, \%options);

        if (defined $vista_feature) {
            push @features_to_store, $vista_feature;
        }

        if (scalar @features_to_store == $batch_size) {
            $adaptors->{'ExternalFeature'}->store(@features_to_store);
            @features_to_store = ();
        }
    }
    $adaptors->{'ExternalFeature'}
             ->store(@features_to_store); #store last batch

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

    my @adaptors_to_get = ('ExternalFeature',
                           'FeatureType',
                           'FeatureSet'
    );

    for my $ad (@adaptors_to_get) {
        $adaptors->{$ad} = $funcgen_db->get_adaptor($ad);
    }

    $adaptors->{'Slice'} = $core_db->get_adaptor('Slice');

    return $adaptors;
}

# Human|chr16:86396481-86397120 | element 1 | positive  | neural tube[12/12] | hindbrain (rhombencephalon)[12/12] | limb[3/12] | cranial nerve[8/12]
sub parse_input_line {
    my ($line, $adaptors, $cache, $options) = @_;

    my @fields = split /\|/, $line;

    shift @fields; # discard first element
    for my $field (@fields) {
        $field =~ s/^\s+|\s+$//g; # trim leading and trailing whitespace
    }

    my $region_string = $fields[0];
    $region_string =~ s/\-/:/;
    my ($coord_system_name, $version);
    if ($region_string =~ /^chr/) {
        $coord_system_name = 'chromosome';
        $region_string =~ s/^chr//g;
    }
    else {
        $coord_system_name = 'scaffold';
    }
    $version   = $options->{assembly};
    my $strand = 1;

    my $name =
        join(':', ($coord_system_name, $version, $region_string, $strand));
    $fields[0] = $name;

    $fields[1] =~ s/\D//g; # discard non-numeric characters in display label

    my $vista_feature =
        create_vista_feature($adaptors, \@fields, $cache);

    return $vista_feature;
}

sub create_vista_feature {
    my ($adaptors, $fields, $cache) = @_;

    my ($name, $display_label, $sign) = @{$fields};

    my $feature_slice = $adaptors->{Slice}->fetch_by_name($name);

    my $slice = $adaptors->{Slice}->fetch_by_region(
        $feature_slice->coord_system->name(),
        $feature_slice->seq_region_name(),
        undef, #start
        undef, #end
        undef, #strand
        $feature_slice->coord_system->version()
    );

    my $feature_type;
    if ($sign eq 'positive') {
        $feature_type = $cache->{positive_feature_type};
    }
    elsif ($sign eq 'negative') {
        $feature_type = $cache->{negative_feature_type};
    }

    my $vista_feature = Bio::EnsEMBL::Funcgen::ExternalFeature->new(
        -feature_set   => $cache->{feature_set},
        -feature_type  => $feature_type,
        -slice         => $slice,
        -start         => $feature_slice->seq_region_start(),
        -end           => $feature_slice->seq_region_end(),
        -display_label => $display_label,
        -strand        => $feature_slice->strand()
    );

    return $vista_feature;
}
