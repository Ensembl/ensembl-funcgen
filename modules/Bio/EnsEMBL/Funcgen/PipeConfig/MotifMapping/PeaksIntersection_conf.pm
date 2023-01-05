=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::PipeConfig::MotifMapping::PeaksIntersection_conf

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

package Bio::EnsEMBL::Funcgen::PipeConfig::MotifMapping::PeaksIntersection_conf;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('Bio::EnsEMBL::Funcgen::PipeConfig::MotifMapping::Base_conf');

sub pipeline_analyses {
    my $self = shift;

    return [
        {
            -logic_name => 'start_PeaksIntersection',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into  => {
                '1->A' => 'matrix_factory_3',
                'A->1' => 'create_MotifFeature_Peak_loading_files',
            },
        },
        {
            -logic_name => 'matrix_factory_3',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -parameters => { inputquery   =>
                             q(select distinct bm.stable_id from binding_matrix bm join binding_matrix_transcription_factor_complex using(binding_matrix_id) join transcription_factor_complex_composition using (transcription_factor_complex_id) join transcription_factor tf using (transcription_factor_id) join peak_calling pc using (feature_type_id)),
                             db_conn      => 'funcgen:#species#',
                             column_names =>
                             [ 'binding_matrix_stable_id' ]
            },
            -flow_into  => {
                2 => { 'Peak_factory' => INPUT_PLUS() },
            }
        },
        {
            -logic_name => 'Peak_factory',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -parameters =>
            { inputquery     => q(select distinct transcription_factor.feature_type_id, peak_calling.epigenome_id from binding_matrix join binding_matrix_transcription_factor_complex using(binding_matrix_id) join transcription_factor_complex_composition using (transcription_factor_complex_id) join transcription_factor using (transcription_factor_id) join peak_calling using (feature_type_id) where binding_matrix.stable_id = '#binding_matrix_stable_id#' ),
              db_conn      => 'funcgen:#species#',
              column_names => [ 'feature_type_id', 'epigenome_id' ]
            },
            -flow_into  => {
                2 => { 'run_Peak_intersection' => INPUT_PLUS() },
            }
        },
        {
            -logic_name => 'run_Peak_intersection',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd => 'bedtools intersect -f 1 -wa -wb -sorted -a '
                    . '#MOODS_bed_dir#/#binding_matrix_stable_id#.sorted.bed'
                    . ' -b '
                    . '#Peaks_dir#/bed/#feature_type_id#__#epigenome_id#.sorted.bed'
                    . ' >#Peaks_intersection_dir#/unfiltered/#binding_matrix_stable_id#___#feature_type_id#__#epigenome_id#.sorted.bed'
            },
            -rc_name    => '2Gb_job',
            -analysis_capacity => 300,
            -batch_size => 10,
            -flow_into  => {
                1 => { 'filter_Peak_intersection' => INPUT_PLUS() },
            }
        },
        {
            -logic_name => 'filter_Peak_intersection',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::'
                . 'FilterPeakIntersection',
            -rc_name    => '1Gb_job',
            -analysis_capacity => 10,
            -batch_size => 200
        },
        {
            -logic_name => 'create_MotifFeature_Peak_loading_files',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::'
                . 'CreateMotifFeaturePeakLoadingFiles',
            -flow_into  => {
                '1' => 'deduplicate',
            },
            -rc_name    => '16Gb_job'
        },
        {
            -logic_name => 'deduplicate',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd =>
                'sort -n #Peaks_intersection_dir#/motif_feature_table.tsv '
                    . '| uniq >#Peaks_intersection_dir#/motif_feature_table.uniq.tsv'
            },
            -rc_name    => '6Gb_job',
            -flow_into  => {
                1 => { 'populate_motif_feature_table' },
            }
        },
        {
            -logic_name => 'populate_motif_feature_table',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
            -parameters => {
                sql     => q(load data local infile ')
                    . '#Peaks_intersection_dir#/motif_feature_table.uniq.tsv'
                    . q(' into table motif_feature fields terminated by '\t' ENCLOSED BY '' lines terminated by '\n' (motif_feature_id, binding_matrix_id,seq_region_id,seq_region_start,seq_region_end,seq_region_strand,score)),
                db_conn => 'funcgen:#species#',
            },
            -flow_into  => {
                1 => { 'populate_motif_feature_peak_table' },
            }
        },
        {
            -logic_name => 'populate_motif_feature_peak_table',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
            -parameters => {
                sql     => q(load data local infile ')
                    . '#Peaks_intersection_dir#/motif_feature_peak_table.tsv'
                    . q(' into table motif_feature_peak fields terminated by '\t' ENCLOSED BY '' lines terminated by '\n' (motif_feature_id, peak_id)),
                db_conn => 'funcgen:#species#',
            },
        },
    ];
}

1;
