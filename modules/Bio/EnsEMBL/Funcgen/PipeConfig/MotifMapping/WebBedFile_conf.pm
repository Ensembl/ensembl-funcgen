=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::PipeConfig::MotifMapping::PeaksIntersection_conf

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

package Bio::EnsEMBL::Funcgen::PipeConfig::MotifMapping::WebBedFile_conf;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('Bio::EnsEMBL::Funcgen::PipeConfig::MotifMapping::Base_conf');

sub pipeline_analyses {
    my $self = shift;

    return [
        {
            -logic_name => 'start_WebBedFile',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into  => {
                '1' => 'create_bed_for_web_track',
            },
        },
        {
            -logic_name => 'create_bed_for_web_track',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::'
                . 'CreateWebBed',
            -rc_name    => '16Gb_job',
            -flow_into  => {
               '1' => 'parse_bed'
            },
        },
        {
            -logic_name => 'parse_bed',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::'
                . 'ParseWebBed',
            -rc_name    => '16Gb_job',
            -flow_into  => {
                '1' => 'bedSort'
            },
        },

        {
            -logic_name => 'bedSort',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => { cmd =>
                             'bedSort  #web_bed_dir#/web_track.parsed.bed #web_bed_dir#/web_track.parsed.sorted.bed'
            },
            -rc_name    => '64Gb_job',
            -flow_into  => {
                1 => 'convert_bed_to_bb',
            }
        },
        {
            -logic_name => 'convert_bed_to_bb',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => { cmd =>
                             'bedToBigBed '
                                . '-tab -as=#auto_sql# -type=bed4+8 '
                                . '-extraIndex=name,binding_matrix_stable_id '
                                . '#web_bed_dir#/web_track.parsed.sorted.bed '
                                . '#chr_sizes# '
                                . '#output_dir#/#species#.#assembly#.motif_features.bb'
            },
            -rc_name    => '64Gb_job',
            -flow_into  => {
                1 => 'update_data_file_table',
            }
        },
        {
            -logic_name => 'update_data_file_table',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
            -parameters => {
                sql     => q(update data_file set path=')
                    . '/funcgen/motif_feature/#release#/#species#.#assembly#.motif_features.bb'
                    . q(' WHERE table_name='motif_feature'),
                db_conn => 'funcgen:#species#',
            },
        },
    ];
}

1;