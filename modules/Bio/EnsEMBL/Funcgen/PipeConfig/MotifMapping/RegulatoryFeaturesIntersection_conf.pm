=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::PipeConfig::MotifMapping::RegulatoryFeaturesIntersection_conf

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

package Bio::EnsEMBL::Funcgen::PipeConfig::MotifMapping::RegulatoryFeaturesIntersection_conf;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('Bio::EnsEMBL::Funcgen::PipeConfig::MotifMapping::Base_conf');

sub pipeline_analyses {
    my $self = shift;

    return [
        {
            -logic_name => 'start_RegulatoryFeaturesIntersection',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into  => {
                '1->A' => 'matrix_factory_4',
                'A->1' => 'move_filtered_files',
            },
        },
        {
            -logic_name => 'matrix_factory_4',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
                -parameters => { inputquery   => q(SELECT stable_id FROM binding_matrix),
                                 db_conn      => 'funcgen:#species#',
                                 column_names => [ 'binding_matrix_stable_id' ]
            },
            -flow_into  => {
                2 => { 'run_RegulatoryFeature_intersection' => INPUT_PLUS() },
            }
        },
        {
            -logic_name        => 'run_RegulatoryFeature_intersection',
            -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters        => {
                cmd => 'bedtools intersect -f 1 -wa -wb -sorted -a '
                    . '#MOODS_bed_dir#/#binding_matrix_stable_id#.sorted.bed'
                    . ' -b '
                    . '#RF_dir#/regulatory_features.sorted.bed'
                    . ' >#RF_intersection_dir#/unfiltered/#binding_matrix_stable_id#___RegulatoryFeatures'
            },
            -rc_name           => '2Gb_job',
            -analysis_capacity => 200,
            -batch_size        => 3,
            -max_retry_count   => 1,
            -flow_into  => {
                1 => { 'junction' => INPUT_PLUS() },
            }
        },
        {
            -logic_name => 'junction',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into  => {
                '1' => [ 'filter_RegulatoryFeature_intersection',
                            'epigenome_factory' ],
            },
        },
        {
            -logic_name      => 'filter_RegulatoryFeature_intersection',
            -module          =>
            'Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::'
                . 'FilterRegulatoryFeatureIntersection',
            -rc_name         => '32Gb_job',
            -analysis_capacity => 300,
            -batch_size      => 3,
            -max_retry_count => 1,
            -flow_into       => {
                1        => 'remove_whitespace',
                MEMLIMIT => 'filter_RegulatoryFeature_intersection_highmem'
            }
        },
        {
            -logic_name      => 'filter_RegulatoryFeature_intersection_highmem',
            -module          =>
            'Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::'
                . 'FilterRegulatoryFeatureIntersection',
            # -module =>    'Bio::EnsEMBL::Hive::RunnableDB::Dummy'
            -rc_name         => '64Gb_job',
            -batch_size      => 3,
            -max_retry_count => 1,
            -flow_into       => {
                1        => 'remove_whitespace',
                # MEMLIMIT => 'filter_RegulatoryFeature_intersection_veryhighmem'
            }
        },
        # {
        #     -logic_name      => 'filter_RegulatoryFeature_intersection_veryhighmem',
        #     -module          =>
        #     'Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::'
        #         . 'FilterRegulatoryFeatureIntersection',
        #     # -module =>    'Bio::EnsEMBL::Hive::RunnableDB::Dummy'
        #     -rc_name         => '64Gb_job',
        #     -batch_size      => 3,
        #     -max_retry_count => 1,
        #     -flow_into       => {
        #         1        => 'remove_whitespace'
        #     }
        # },
        {
            -logic_name => 'remove_whitespace',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },

        {
            -logic_name => 'epigenome_factory',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            # TODO join to regulatory_build_epigenome table
            -parameters => { inputquery   => q(select distinct peak_calling.epigenome_id from binding_matrix join binding_matrix_transcription_factor_complex using(binding_matrix_id) join transcription_factor_complex_composition using (transcription_factor_complex_id) join transcription_factor using (transcription_factor_id) join peak_calling using (feature_type_id) where binding_matrix.stable_id ='#binding_matrix_stable_id#'),
                             db_conn      => 'funcgen:#species#',
                             column_names => [ 'epigenome_id' ]
            },
            -flow_into  => {
                2 => { 'filter_RegulatoryFeature_intersection' => INPUT_PLUS() },
            }
        },
        {
            -logic_name        => 'move_filtered_files',
            -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters        => {
                cmd => 'mv #RF_intersection_dir#/unfiltered/*.filtered '
                    . '#RF_intersection_dir#/filtered/'
            },
            -rc_name           => '2Gb_job',
            -max_retry_count   => 1,
            -flow_into  => {
                1 => { 'create_MF_RF_loading_files' => INPUT_PLUS() },
            }
        },

        {
            -logic_name => 'create_MF_RF_loading_files',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::'
                . 'CreateMotifFeatureRegulatoryFeatureLoadingFiles',
            -flow_into  => {
                '1' => 'deduplicate_2',
            },
            -rc_name    => '16Gb_job'
        },
        {
            -logic_name => 'deduplicate_2',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd =>
                'sort -n #RF_intersection_dir#/motif_feature_table.tsv '
                    . '| uniq >#RF_intersection_dir#/motif_feature_table.uniq.tsv'
            },
            -rc_name    => '16Gb_job',
            -flow_into  => {
                1 => { 'populate_motif_feature_table_2' },
            }
        },
        {
            -logic_name => 'populate_motif_feature_table_2',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
            -parameters => {
                sql     => q(load data local infile ')
                    . '#RF_intersection_dir#/motif_feature_table.uniq.tsv'
                    . q(' into table motif_feature fields terminated by '\t' ENCLOSED BY '' lines terminated by '\n' (motif_feature_id, binding_matrix_id,seq_region_id,seq_region_start,seq_region_end,seq_region_strand,score)),
                db_conn => 'funcgen:#species#',
            },
            -flow_into  => {
                1 => { 'populate_motif_feature_regulatory_feature_table' },
            }
        },
        {
            -logic_name => 'populate_motif_feature_regulatory_feature_table',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
            -parameters => {
                sql     => q(load data local infile ')
                    . '#RF_intersection_dir#/motif_feature_regulatory_feature_table.tsv'
                    . q(' into table motif_feature_regulatory_feature fields terminated by '\t' ENCLOSED BY '' lines terminated by '\n' (motif_feature_id, regulatory_feature_id, epigenome_id, has_matching_Peak)),
                db_conn => 'funcgen:#species#',
            },
        },
    ];
}

1;
