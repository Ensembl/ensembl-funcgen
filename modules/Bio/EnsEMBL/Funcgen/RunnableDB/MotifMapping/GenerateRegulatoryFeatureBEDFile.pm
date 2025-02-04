=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::GenerateRegulatoryFeatureBEDFile

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

package Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::GenerateRegulatoryFeatureBEDFile;

use strict;
use autodie;
use base ('Bio::EnsEMBL::Hive::Process');

sub run {
    my $self = shift;

    my $species = $self->param('species');

    my $funcgen_adaptor
        = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');

    my $RF_dir = $self->param_required ('RF_dir');

    my $regulatory_feature_adaptor =
        $funcgen_adaptor->get_RegulatoryFeatureAdaptor;

    #   p $regulatory_feature_adaptor;
    my $regulatory_build_adaptor = $funcgen_adaptor->get_RegulatoryBuildAdaptor;
    my $regulatory_build =
        $regulatory_build_adaptor->fetch_current_regulatory_build;

    my $iterator =
        $regulatory_feature_adaptor->fetch_Iterator_by_RegulatoryBuild(
            $regulatory_build);

    my $regulatory_features_file = $RF_dir . '/regulatory_features.bed';

    open my $regulatory_features_fh, '>', $regulatory_features_file;

    REGULATORY_FEATURE:
    while ( $iterator->has_next ) {

        my $current_regulatory_feature = $iterator->next;

        my $bed_line = join "\t",
                            (
                                $current_regulatory_feature->slice->seq_region_name,
                                $current_regulatory_feature->bound_seq_region_start,
                                $current_regulatory_feature->bound_seq_region_end,
                                $current_regulatory_feature->dbID,
                            );

        $regulatory_features_fh->print($bed_line);
        $regulatory_features_fh->print("\n");
    }
    $regulatory_features_fh->close;

}

1;
