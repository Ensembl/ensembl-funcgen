=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::FilterPeakIntersection

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

package Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::ParseWebBed;

use strict;
use autodie;
use base ('Bio::EnsEMBL::Hive::Process');
use Bio::EnsEMBL::Utils::SqlHelper;

sub run {
    my $self = shift;

    my $species = $self->param_required('species');
    my $assembly = $self->param_required('assembly');
    my $web_bed_dir = $self->param_required('web_bed_dir');

    my $funcgen_adaptor
        = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
    my $core_adaptor =
        Bio::EnsEMBL::Registry->get_DBAdaptor( $species , 'core' );

    my $bma = $funcgen_adaptor->get_adaptor('BindingMatrix');
    my $mfa = $funcgen_adaptor->get_adaptor('MotifFeature');

    my $helper = Bio::EnsEMBL::Utils::SqlHelper->new(
    -DB_CONNECTION => $funcgen_adaptor->dbc
    );

    my $rows = $helper->execute_simple(
    -SQL => 'SELECT MAX(motif_feature_id) FROM motif_feature_peak'
    );

    my $max_experimentally_verified_motif_feature_id = $rows->[0];

    my %tfc_name_cache;

    open my $unparsed_bed_fh, '<', $web_bed_dir . '/web_track.bed';
    open my $bed_hf, '>', $web_bed_dir . '/web_track.parsed.bed';

    while ( readline $unparsed_bed_fh ) {
        chomp;

        my ( $chrom, $start, $end, $mf_stable_id, $score,
            $strand, $thStart, $thEnd, $rgb, $bm_stable_id, $motif_feature_id )
            = split /\t/;

        if ( !exists $tfc_name_cache{$bm_stable_id} ) {
            my $bm        = $bma->fetch_by_stable_id($bm_stable_id);
            my $tfc_names = $bm->get_TranscriptionFactorComplex_names();
            my $tfc_label = join ',', @{$tfc_names};
            $tfc_name_cache{$bm_stable_id} = $tfc_label;
        }

        $rgb = '153,153,153';
        my $epigenome_names;

        if ( $motif_feature_id <= $max_experimentally_verified_motif_feature_id ) {
            $rgb = '0,0,0';
            my $mf             = $mfa->fetch_by_dbID($motif_feature_id);
            my $epigenome_list = $mf->get_all_Epigenomes_with_experimental_evidence;
            my @epigenome_names_list;
            for my $epigenome ( @{$epigenome_list} ) {
                push @epigenome_names_list, $epigenome->short_name;
            }
            $epigenome_names = join ',', @epigenome_names_list;
        }
        else {
            $epigenome_names = 'N/A';
        }

        $bed_hf->print ( $chrom . "\t"
                . $start . "\t"
                . $end . "\t"
                . $mf_stable_id . "\t"
                . $score . "\t"
                . $strand . "\t"
                . $thStart . "\t"
                . $thEnd . "\t"
                . $rgb . "\t"
                . $bm_stable_id . "\t"
                . $tfc_name_cache{$bm_stable_id} . "\t"
                . $epigenome_names . "\n");
    }
    close $unparsed_bed_fh;
}

1;
