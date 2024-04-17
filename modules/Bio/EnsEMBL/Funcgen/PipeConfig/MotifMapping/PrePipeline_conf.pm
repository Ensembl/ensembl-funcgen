=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::PipeConfig::MotifMapping::PrePipeline_conf

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

package Bio::EnsEMBL::Funcgen::PipeConfig::MotifMapping::PrePipeline_conf;

use base ('Bio::EnsEMBL::Funcgen::PipeConfig::MotifMapping::Base_conf');

sub pipeline_analyses {
    my $self = shift;

    return [
        {
            -logic_name => 'start_pre_pipeline_checks_tasks',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into  => {
                MAIN => 'TODO_pre_pipeline_checks',
            },
        },
        {
            -logic_name => 'TODO_pre_pipeline_checks',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into  => {
                MAIN => 'make_dirs',
            },
        },
        {
            -logic_name => 'make_dirs',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd => 'mkdir -p'
                    . ' #matrices_dir#'
                    . ' #MOODS_raw_dir#'
                    . ' #MOODS_bed_dir#'
                    . ' #Peaks_dir#'
                    . ' #RF_dir#'
                    . ' #Peaks_intersection_dir#/filtered'
                    . ' #Peaks_intersection_dir#/unfiltered'
                    . ' #RF_intersection_dir#/filtered'
                    . ' #RF_intersection_dir#/unfiltered'
                    . ' #stable_id_mapping_dir#'
                    . ' #temp_dir#'
                    . ' #output_dir#'
                    . ' #web_bed_dir#'
                    . ' #ftp_dumps_dir#'
                    . ' #ftp_dumps_dir#/chromosomal'
                    . ' #ftp_dumps_dir#/non_chromosomal'
            },
        }
    ];
}

1;