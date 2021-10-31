=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::FilterPeakIntersection

=head1 DESCRIPTION



=head1 LICENSE

    Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
    Copyright [2016-2021] EMBL-European Bioinformatics Institute

    Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

         http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software distributed under the License
    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and limitations under the License.

=head1 CONTACT

    ensembl-funcgen@ebi.ac.uk

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::ExportMotifFeaturesToBed;

use strict;
use base ('Bio::EnsEMBL::Hive::Process');

sub run {
    my $self = shift;

    my $species = $self->param_required('species');
    my $stable_id_mapping_dir = $self->param_required('stable_id_mapping_dir');

    my $funcgen_adaptor
        = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
    my $previous_funcgen_adaptor =
        Bio::EnsEMBL::Registry->get_DBAdaptor( $species . '_previous_version',
                                               'funcgen' );

    export($funcgen_adaptor, $stable_id_mapping_dir, 'current');
    export($previous_funcgen_adaptor, $stable_id_mapping_dir, 'previous');

}

sub export {
    my ($dba, $stable_id_mapping_dir, $suffix) = @_;

    my $motif_feature_adaptor = $dba->get_adaptor('MotifFeature');

    my $outfile = $stable_id_mapping_dir . '/' . 'motif_features_'
        . $suffix . '.bed';

    open my $out_fh, '>', $outfile;
    $motif_feature_adaptor->_bulk_export_to_bed($out_fh);
    close $out_fh;
}


1;
