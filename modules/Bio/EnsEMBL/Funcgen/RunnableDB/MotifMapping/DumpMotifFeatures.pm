=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::FilterPeakIntersection

=head1 DESCRIPTION



=head1 LICENSE

    Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
    Copyright [2016-2023] EMBL-European Bioinformatics Institute

    Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

         http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software distributed under the License
    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and limitations under the License.

=head1 CONTACT

    ensembl-funcgen@ebi.ac.uk

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::DumpMotifFeatures;

use strict;
use base ('Bio::EnsEMBL::Hive::Process');
use autodie;
use IO::File;
use Bio::EnsEMBL::Utils::IO::GFFSerializer;
use Bio::EnsEMBL::Utils::SqlHelper;

sub run {
    my $self = shift;

    my $species = $self->param_required('species');

    my $funcgen_adaptor
        = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
    my $ontology_term_adaptor =
        Bio::EnsEMBL::Registry->get_adaptor('Multi', 'Ontology', 'OntologyTerm');

    my $ftp_dumps_dir = $self->param_required('ftp_dumps_dir');

    my $output_file = $ftp_dumps_dir . '/motif_features.gff';

    my $output_fh = IO::File->new(">$output_file");

    my $serializer = Bio::EnsEMBL::Utils::IO::GFFSerializer->new(
        $ontology_term_adaptor,
        $output_fh
    );

    my $motif_feature_adaptor =
        $funcgen_adaptor->get_adaptor('MotifFeature');

    my $core_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'Core' );
    $core_adaptor->dbc->disconnect_when_inactive(0);
    $ontology_term_adaptor->db->dbc->disconnect_when_inactive(0);
    $funcgen_adaptor->dbc->disconnect_when_inactive(0);

    my $batch_size = 10_000;
    my @fetch_list;

    my $mf_id_cnt = 1;

    my $helper = Bio::EnsEMBL::Utils::SqlHelper->new(
        -DB_CONNECTION => $funcgen_adaptor->dbc
    );

    my $total_number_of_motif_features = $helper->execute_simple(
        -SQL      => 'select count(motif_feature_id) from motif_feature',
    )->[0];


    while ($mf_id_cnt <= $total_number_of_motif_features) {
        push @fetch_list, $mf_id_cnt;

        if ($mf_id_cnt % $batch_size == 0 || $mf_id_cnt == $total_number_of_motif_features) {
            my $motif_features =
                $motif_feature_adaptor->fetch_all_by_dbID_list(\@fetch_list);
            for my $motif_feature (@{$motif_features}) {
                $serializer->print_feature($motif_feature);
            }
            @fetch_list = ();
        }
        $mf_id_cnt++;
    }

}

1;
