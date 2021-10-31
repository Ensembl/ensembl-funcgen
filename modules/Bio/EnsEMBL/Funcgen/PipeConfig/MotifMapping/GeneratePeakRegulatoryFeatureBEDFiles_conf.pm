=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::PipeConfig::MotifMapping::GeneratePeakRegulatoryFeatureFiles_conf

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

package Bio::EnsEMBL::Funcgen::PipeConfig::MotifMapping::GeneratePeakRegulatoryFeatureBEDFiles_conf;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('Bio::EnsEMBL::Funcgen::PipeConfig::MotifMapping::Base_conf');

sub pipeline_analyses {
    my $self = shift;

    return [
        { #
            -logic_name => 'generate_Peak_RF_BED_Files',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into  => {
                '1' => [ 'Peak_BED_file_factory',
                         'generate_RF_BED_file' ],
            },
        },
        {
            -logic_name => 'generate_RF_BED_file',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping'
                . '::GenerateRegulatoryFeatureBEDFile',
            -parameters => {},
            -flow_into  => {
                1 => { 'sort_RF_BED_file' },
            }
        },
        {
            -logic_name => 'sort_RF_BED_file',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd =>
                'sort -k1,1 -k2,2n -k3,3n #RF_dir#/regulatory_features.bed >#RF_dir#/regulatory_features.sorted.bed'
            },
            -rc_name    => '16Gb_job',
        },

        {
            -logic_name => 'Peak_BED_file_factory',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -parameters =>
            { inputquery     => q(select distinct transcription_factor.feature_type_id, peak_calling.epigenome_id from binding_matrix join binding_matrix_transcription_factor_complex using(binding_matrix_id) join transcription_factor_complex_composition using (transcription_factor_complex_id) join transcription_factor using (transcription_factor_id) join peak_calling using (feature_type_id)),
              db_conn      => 'funcgen:#species#',
              column_names => [ 'feature_type_id', 'epigenome_id' ]
            },
            -flow_into  => {
                2 => { 'generate_Peak_BED_files' => INPUT_PLUS() },
            }
        },
        {
            -logic_name => 'generate_Peak_BED_files',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping'
                . '::GeneratePeakBEDFiles',
            -parameters => {},
            # -rc_name    => '16Gb_job',
            -flow_into  => {
                1 => { 'sort_Peak_BED_files' => INPUT_PLUS() },
            }
        },
        {
            -logic_name        => 'sort_Peak_BED_files',
            -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters        => {
                cmd =>
                'sort -k1,1 -k2,2n -k3,3n #Peaks_dir#/bed/#feature_type_id#__#epigenome_id#.bed >#Peaks_dir#/bed/#feature_type_id#__#epigenome_id#.sorted.bed'
            },
            # -rc_name           => '16Gb_job',
            -analysis_capacity => 300
        }
    ];
}

1;
