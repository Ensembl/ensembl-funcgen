=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::PipeConfig::MotifMapping::BindingMatrixTranscriptionFactorRegistration_conf

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

package Bio::EnsEMBL::Funcgen::PipeConfig::MotifMapping::BindingMatrixTranscriptionFactorRegistration_conf;

use base ('Bio::EnsEMBL::Funcgen::PipeConfig::MotifMapping::Base_conf');

sub pipeline_analyses {
    my $self = shift;

    return [
        {
            -logic_name => 'start_registration',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into  => {
                MAIN => 'truncate_tables',
            },
        },
        {
            -logic_name => 'truncate_tables',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
            -parameters => {
                sql     => [
                    "TRUNCATE binding_matrix;",
                    "TRUNCATE binding_matrix_frequencies;",
                    "TRUNCATE binding_matrix_transcription_factor_complex;",
                    "TRUNCATE transcription_factor_complex_composition;",
                    "TRUNCATE transcription_factor_complex;",
                    "TRUNCATE transcription_factor;",
                    "TRUNCATE motif_feature;",
                    "TRUNCATE motif_feature_peak;",
                    "TRUNCATE motif_feature_regulatory_feature;",
                ],
                db_conn => 'funcgen:#species#',
            },
            -flow_into  => {
                MAIN => 'register_binding_matrices_transcription_factors',
            },
        },
        {
            -logic_name => 'register_binding_matrices_transcription_factors',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::'
                . 'RegisterBindingMatrixTranscriptionFactor',

        }
    ];
}

1;