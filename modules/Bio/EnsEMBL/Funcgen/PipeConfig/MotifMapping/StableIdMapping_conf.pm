=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::PipeConfig::MotifMapping::PeaksIntersection_conf

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

package Bio::EnsEMBL::Funcgen::PipeConfig::MotifMapping::StableIdMapping_conf;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('Bio::EnsEMBL::Funcgen::PipeConfig::MotifMapping::Base_conf');

sub pipeline_analyses {
    my $self = shift;

    return [
        {
            -logic_name => 'start_StableIdMapping',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into  => {
                '1' => 'export_motif_features_to_bed',
            },
        },
        {
            -logic_name => 'export_motif_features_to_bed',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::'
                . 'ExportMotifFeaturesToBed',
            -rc_name    => '64Gb_job',
            -flow_into  => {
               '1->A' => ['sort_current_bed',
                         'sort_previous_bed'],
               'A->1' => 'intersect'
            },
        },
        {
            -logic_name => 'sort_current_bed',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd => 'sort --parallel=16 -T #temp_dir# --buffer-size=31G '
                . '-k1,1 -k2,2n -k3,3n '
                . '<#stable_id_mapping_dir#/motif_features_current.bed '
                . '>#stable_id_mapping_dir#/motif_features_current.sorted.bed'
            },
            -rc_name    => '16c32Gb_job',
        },
        {
            -logic_name => 'sort_previous_bed',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd => 'sort --parallel=16 -T #temp_dir# --buffer-size=31G '
                    . '-k1,1 -k2,2n -k3,3n '
                    . '<#stable_id_mapping_dir#/motif_features_previous.bed '
                    . '>#stable_id_mapping_dir#/motif_features_previous.sorted.bed'
            },
            -rc_name    => '16c32Gb_job',
        },
        {
            -logic_name => 'intersect',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd =>
                'bedtools intersect -loj -sorted -s -f 1 -F 1 -wa -wb '
                . '-a #stable_id_mapping_dir#/motif_features_current.sorted.bed '
                . '-b #stable_id_mapping_dir#/motif_features_previous.sorted.bed '
                . '>#stable_id_mapping_dir#/intersection.out'
            },
            -rc_name    => '32Gb_job',
            -flow_into  => {
                '1' => 'sort_intersection',
            },
        },
        {
            -logic_name => 'sort_intersection',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd => 'sort --parallel=16 -T #temp_dir# --buffer-size=31G '
                . '-k4,4n '
                . '<#stable_id_mapping_dir#/intersection.out '
                . '>#stable_id_mapping_dir#/intersection.sorted.out'
            },
            -rc_name    => '16c32Gb_job',
            -flow_into => {
                '1' => 'assign_stable_ids'
            }
        },
        {
            -logic_name => 'assign_stable_ids',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::'
                . 'AssignStableIDs',
            -rc_name    => '32Gb_job',
            -flow_into  => {
                '1' => 'create_temp_table'
            },
        },
        {
            -logic_name => 'create_temp_table',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
            -parameters => {
                sql     => q(CREATE TABLE motif_feature_temp LIKE motif_feature),
                db_conn => 'funcgen:#species#',
            },
            -flow_into  => {
                '1' => 'alter_temp_table'
            },
        },
        {
            -logic_name => 'alter_temp_table',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
            -parameters => {
                sql     => 'ALTER TABLE motif_feature_temp '
                            . 'DROP COLUMN score, '
                            . 'DROP COLUMN seq_region_strand, '
                            . 'DROP COLUMN seq_region_end, '
                            . 'DROP COLUMN seq_region_start, '
                            . 'DROP COLUMN seq_region_id, '
                            . 'DROP COLUMN binding_matrix_id, '
                            . 'DROP INDEX binding_matrix_idx, '
                            . 'DROP INDEX seq_region_idx, '
                            . 'DROP INDEX unique_idx',
                db_conn => 'funcgen:#species#',
            },
            -flow_into  => {
                '1' => 'populate_temp_table'
            },
        },
        {
            -logic_name => 'populate_temp_table',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
            -parameters => {
                sql     => q(load data local infile ')
                    . '#stable_id_mapping_dir#/new_stableIDs.out'
                    . q(' into table motif_feature_temp fields terminated by '\t' ENCLOSED BY '' lines terminated by '\n' (motif_feature_id, stable_id)),
                db_conn => 'funcgen:#species#',
            },
            -flow_into  => {
                '1' => 'update_stableIDs'
            },
        },
        {
            -logic_name => 'update_stableIDs',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
            -parameters => {
                sql     => 'update motif_feature a join motif_feature_temp b using(motif_feature_id) set a.stable_id = b.stable_id',
                db_conn => 'funcgen:#species#',
            },
            -flow_into  => {
                '1' => 'drop_temp_table'
            },
        },
        {
            -logic_name => 'drop_temp_table',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
            -parameters => {
                sql     => 'DROP TABLE motif_feature_temp',
                db_conn => 'funcgen:#species#',
            },
        },
    ];
}

1;