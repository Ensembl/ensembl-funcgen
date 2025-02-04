=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::FilterRegulatoryFeatureIntersection

=head1 DESCRIPTION



=head1 LICENSE

    Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
    Copyright [2016-2025] EMBL-European Bioinformatics Institute

    Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

         http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software distributed under the License
    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and limitations under the License.

=head1 CONTACT

    ensembl-funcgen@ebi.ac.uk

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::FilterRegulatoryFeatureIntersection;

use strict;
use autodie;
use base ('Bio::EnsEMBL::Hive::Process');

sub run {
    my $self = shift;

    my $species             = $self->param_required('species');
    my $RF_intersection_dir = $self->param_required('RF_intersection_dir');
    my $binding_matrix_stable_id              = $self->param_required('binding_matrix_stable_id');
    my $epigenome_id      = $self->param('epigenome_id');

    my $funcgen_adaptor =
        Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
    my $core_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'core');

    # load all experimentally verified Motif Features for a given Binding Matrix and Epigenome
    my $experimentally_verified_Motif_Features;
    my $epigenome;
    if ($epigenome_id) {
        $experimentally_verified_Motif_Features =
            $self->get_all_experimentally_verified_Motif_Features();
        my $epigenome_adaptor = $funcgen_adaptor->get_adaptor('Epigenome');

        $epigenome = $epigenome_adaptor->fetch_by_dbID($epigenome_id);
    }

    # create caches for minimizing db connections
    my $chr2seqRegionID; # cache for converting chromosome names to region IDs
    my $seqRegionID2chr; # cache for converting region IDs to chromosome names

    my $slice_adaptor = $core_adaptor->get_adaptor('Slice');
    my $slices        = $slice_adaptor->fetch_all('toplevel');
    for my $slice (@{$slices}) {
        $chr2seqRegionID->{ $slice->seq_region_name }   =
            $slice->get_seq_region_id;
        $seqRegionID2chr->{ $slice->get_seq_region_id } = $slice->seq_region_name;
    }

    my $regulatory_feature_adaptor =
        $funcgen_adaptor->get_adaptor('RegulatoryFeature');
    my $binding_matrix_adaptor = $funcgen_adaptor->get_adaptor('BindingMatrix');
    my $binding_matrix =
        $binding_matrix_adaptor->fetch_by_stable_id($binding_matrix_stable_id);
    my $binding_matrix_id = $binding_matrix->dbID;

    my $filename = $RF_intersection_dir
        . '/unfiltered/'
        . $binding_matrix_stable_id
        . '___RegulatoryFeatures';

    # load motif feature / regulatory feature overlaps in a hash, grouped by regulatory_feature_id
    my $overlaps = group_overlaps_by_rfID($filename);

    filter($overlaps, $regulatory_feature_adaptor, $binding_matrix_id, $chr2seqRegionID, $filename, $epigenome, $experimentally_verified_Motif_Features);
}

sub filter {
    my ($overlaps, $regulatory_feature_adaptor, $binding_matrix_id, $chr2seqRegionID, $filename, $epigenome, $experimentally_verified_Motif_Features) = @_;

    for my $rf_id ( keys %{$overlaps} ) {
        my $regulatory_feature =
            $regulatory_feature_adaptor->fetch_by_dbID($rf_id);

        if ($epigenome) {

            my $rf_activity_obj =
                $regulatory_feature->regulatory_activity_for_epigenome($epigenome);

            my $rf_activity;
            if ($rf_activity_obj){
                $rf_activity=$rf_activity_obj->activity;
            }

            # skip inactive regulatory features
            if (  !$rf_activity
                || $rf_activity eq 'NA'
                || $rf_activity eq 'INACTIVE' )
            {
                next;
            }

            for
            my $overlap ( @{ sort_overlaps_by_score( $overlaps->{$rf_id} ) } )
            {
                my $is_experimentally_verified =
                    is_experimentally_verified($overlap, $binding_matrix_id, $chr2seqRegionID, $experimentally_verified_Motif_Features);
                if ($is_experimentally_verified) {
                    print_overlap_line( $overlap, $rf_id, $filename, $epigenome );
                    last;
                }
            }

        }
        else {
            my $highest_scoring_overlap =
                sort_overlaps_by_score( $overlaps->{$rf_id} )->[0];
            print_overlap_line( $highest_scoring_overlap, $rf_id, $filename );
        }

    }

}

sub get_all_experimentally_verified_Motif_Features {
    my $self                  = shift;
    my $Peak_intersection_dir = $self->param_required('Peaks_intersection_dir');

    open my $mf_table_fh, '<',
         $Peak_intersection_dir . '/motif_feature_table.uniq.tsv';

    my %experimentally_verified_Motif_Features;

    while (readline $mf_table_fh) {
        my $line = $_;
        chomp $line;

        my @fields            = split /\t/, $line;
        my $binding_matrix_id = $fields[1];
        my $seq_region_id     = $fields[2];
        my $seq_region_start  = $fields[3];
        my $seq_region_strand = $fields[5];

        my $unique_key =
            $binding_matrix_id . '_'
                . $seq_region_id . '_'
                . $seq_region_start . '_'
                . $seq_region_strand;

        $experimentally_verified_Motif_Features{$unique_key} = 1;

    }
    return \%experimentally_verified_Motif_Features;
}

sub group_overlaps_by_rfID {
    my ($filename) = @_;
    my $overlaps = {};

    open my $unfiltered_overlap_fh, '<', $filename;

    while ( readline $unfiltered_overlap_fh ) {
        my $line = $_;
        chomp $line;
        my @fields = split /\t/, $line;
        my $rf_id = $fields[9];
        my ( $mf_score, $seq_region_name ) = split "_", $fields[3];
        my $seq_region_strand = $fields[5];
        if ( $seq_region_strand eq '+' ) {
            $seq_region_strand = 1;
        }
        elsif ( $seq_region_strand eq '-' ) {
            $seq_region_strand = -1;
        }
        else {
            die 'Invalid seq_region_strand: ' . $seq_region_strand;
        }

        $overlaps->{$rf_id} //= [];
        push @{ $overlaps->{$rf_id} },
             {
                 'chr'             => $fields[0],
                 'mf_start'        => $fields[1],
                 'mf_end'          => $fields[2],
                 'mf_score'        => $mf_score,
                 'seq_region_name' => $seq_region_name,
                 'strand'          => $seq_region_strand,
                 'rf_start'        => $fields[7],
                 'rf_end'          => $fields[8],
             };

    }
    return $overlaps;
}

sub is_experimentally_verified {
    my ($overlap, $binding_matrix_id, $chr2seqRegionID, $experimentally_verified_Motif_Features) = @_;

    my $unique_key =
        $binding_matrix_id . '_'
            . $chr2seqRegionID->{ $overlap->{chr} } . '_'
            . $overlap->{mf_start} . '_'
            . $overlap->{strand};


    if ( exists $experimentally_verified_Motif_Features->{$unique_key} ) {
        return 1;
    }
    else {
        return 0;
    }
}

sub print_overlap_line {
    my ( $overlap, $rf_id, $filename, $epigenome ) = @_;

    my $filtered_file = $filename;

    if ($epigenome) {
        $filtered_file .= '_' . $epigenome->short_name . '.filtered';
    }
    else {
        $filtered_file .= '_no_epigenome.filtered';
    }

    open my $filtered_fh, '>>', $filtered_file;

    my $line =
        $overlap->{chr} . "\t"
            . $overlap->{mf_start} . "\t"
            . $overlap->{mf_end} . "\t"
            . $overlap->{mf_score} . "\t"
            . $overlap->{strand} . "\t"
            . $overlap->{seq_region_name} . "\t"
            . $overlap->{rf_start} . "\t"
            . $overlap->{rf_end} . "\t"
            . $rf_id;

    if ($epigenome) {
        $line .= "\t" . $epigenome->dbID;
    }

    $line .= "\n";

    print $filtered_fh $line;

}

sub sort_overlaps_by_score {
    my ($unsorted_overlaps) = @_;

    my @sorted_overlaps =
        sort { $b->{'mf_score'} <=> $a->{'mf_score'} } @{$unsorted_overlaps};

    return \@sorted_overlaps;
}

1;