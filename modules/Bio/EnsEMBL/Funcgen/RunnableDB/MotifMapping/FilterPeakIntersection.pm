=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::FilterPeakIntersection

=head1 DESCRIPTION



=head1 LICENSE

    Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
    Copyright [2016-2024] EMBL-European Bioinformatics Institute

    Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

         http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software distributed under the License
    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and limitations under the License.

=head1 CONTACT

    ensembl-funcgen@ebi.ac.uk

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::FilterPeakIntersection;

use strict;
use base ('Bio::EnsEMBL::Hive::Process');

sub run {
    my $self = shift;

    # my $species = $self->param('species');
    #
    # my $funcgen_adaptor
    #     = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');

    my $Peaks_intersection_dir     = $self->param_required('Peaks_intersection_dir');
    my $binding_matrix_stable_id        = $self->param_required('binding_matrix_stable_id');
    my $feature_type_id = $self->param_required('feature_type_id');
    my $epigenome_id = $self->param_required('epigenome_id');

    open my $unfiltered_intersection_fh, '<',
         $Peaks_intersection_dir
             . '/unfiltered/'
             . $binding_matrix_stable_id
             . '___'
             . $feature_type_id
             . '__'
             . $epigenome_id
             . '.sorted.bed';

    my %filtered_overlaps;
    while (readline $unfiltered_intersection_fh) {
        my @fields        = split /\t/;
        my ($motif_score) = split "_", $fields[3];
        my $peak_id       = $fields[9];
        if (!exists $filtered_overlaps{$peak_id}) {
            $filtered_overlaps{$peak_id} = \@fields;
        }
        else {
            my ($max_motif_score) =
                split "_", $filtered_overlaps{$peak_id}->[3];
            if ($motif_score > $max_motif_score) {
                $filtered_overlaps{$peak_id} = \@fields;
            }
        }
    }

    close $unfiltered_intersection_fh;

    open my $filtered_intersection_fh, '>',
         $Peaks_intersection_dir
             . '/filtered/'
             . $binding_matrix_stable_id
             . '___'
             . $feature_type_id
             . '__'
             . $epigenome_id
             . '.filtered.bed';

    for my $value (values %filtered_overlaps) {
        print $filtered_intersection_fh join("\t", @{$value});
    }

    close $filtered_intersection_fh;
}

1;
