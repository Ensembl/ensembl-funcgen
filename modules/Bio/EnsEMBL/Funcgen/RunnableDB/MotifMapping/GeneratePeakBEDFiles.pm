=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::GeneratePeakBEDFiles

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

package Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::GeneratePeakBEDFiles;

use strict;
use base ('Bio::EnsEMBL::Hive::Process');

sub run {
    my $self = shift;

    my $species = $self->param('species');

    my $funcgen_adaptor
        = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');

    my $Peaks_dir = $self->param_required('Peaks_dir');

    my $epigenome_id = $self->param_required('epigenome_id');
    my $feature_type_id = $self->param_required('feature_type_id');
    my $peak_adaptor = $funcgen_adaptor->get_adaptor('Peak');
    my $peak_calling_adaptor = $funcgen_adaptor->get_adaptor('PeakCalling');

    my $output_dir = $Peaks_dir . '/bed';

    my $cmd = 'mkdir -p ' . $output_dir;
    system($cmd);


    my $feature_type_adaptor = $funcgen_adaptor->get_adaptor('FeatureType');
    my $epigenome_adaptor = $funcgen_adaptor->get_adaptor('Epigenome');
    my $epigenome    = $epigenome_adaptor->fetch_by_dbID($epigenome_id);
    my $feature_type = $feature_type_adaptor->fetch_by_dbID($feature_type_id);
    my $peak_calling = shift @{$peak_calling_adaptor->fetch_all_by_Epigenome_FeatureType($epigenome, $feature_type)};

    my $bed_name = $feature_type_id . '__' . $epigenome_id . '.bed';

    open my $bed_fh, '>', $output_dir . '/' . $bed_name;
    $peak_adaptor->_bulk_export_to_bed_by_PeakCalling($peak_calling,
                                                      $bed_fh);
    close $bed_fh;
}

1;
