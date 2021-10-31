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

package Bio::EnsEMBL::Funcgen::PipeConfig::MotifMapping::FtpDumps_conf;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('Bio::EnsEMBL::Funcgen::PipeConfig::MotifMapping::Base_conf');

sub pipeline_analyses {
    my $self = shift;

    return [
        {
            -logic_name => 'start_FtpDumps',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into  => {
                '1' => 'dump_motif_features_to_gff',
            },
        },
        {
            -logic_name => 'dump_motif_features_to_gff',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::'
                . 'DumpMotifFeatures',
            -rc_name    => '16Gb_job',
            -flow_into  => {
               '1' => 'sort_dumps'
            },
        },
        {
            -logic_name => 'sort_dumps',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd => 'sort --parallel=16 -T #temp_dir# --buffer-size=31G '
                    . '-k1,1 -k4,4n -k5,5n '
                    . '<#ftp_dumps_dir#/motif_features.gff '
                    . '>#ftp_dumps_dir#/motif_features.sorted.gff'
            },
            -rc_name    => '16c32Gb_job',
            -flow_into => {
                '1' => 'split_to_chrs'
            }
        },
        #awk '{print>$1".gff"}' ../Homo_sapiens.GRCh38.motif_features.sorted.gff
        {
            -logic_name => 'split_to_chrs',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd => q(awk '{print>"#ftp_dumps_dir#/#species#.#assembly#."$1".motif_features.gff"}' )
                    . '#ftp_dumps_dir#/motif_features.sorted.gff'
            },
            -rc_name    => '16Gb_job',
            -flow_into => {
                '1' => 'move_files'
            }
        },
        {
            -logic_name => 'move_files',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd => q(mv #ftp_dumps_dir#/#species#.#assembly#.[0-9]* #ftp_dumps_dir#/chromosomal && )
                        . q(mv #ftp_dumps_dir#/#species#.#assembly#.MT* #ftp_dumps_dir#/chromosomal &&)
                        . q(mv #ftp_dumps_dir#/#species#.#assembly#.X.* #ftp_dumps_dir#/chromosomal &&)
                        . q(mv #ftp_dumps_dir#/#species#.#assembly#.Y.* #ftp_dumps_dir#/chromosomal &&)
                        . q(mv #ftp_dumps_dir#/#species#.#assembly#* #ftp_dumps_dir#/non_chromosomal)

            },
            -flow_into => {
                '1' => [ 'bgzip',
                         'gzip_chromosomal',
                         'gzip_non_chromosomal' ],
            }
        },
        {
            -logic_name => 'bgzip',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd => q(bgzip <#ftp_dumps_dir#/motif_features.sorted.gff >#output_dir#/#species#.#assembly#.motif_features.gff.gz )
            },
            -flow_into => {
                '1' => 'tabix'
            },
            -rc_name    => '32Gb_job'
        },
        {
            -logic_name => 'gzip_chromosomal',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd => q(gzip #ftp_dumps_dir#/chromosomal/*gff && )
                        . q(mv #ftp_dumps_dir#/chromosomal/*gz #output_dir#  )
            },
            -rc_name    => '32Gb_job'
        },
        {
            -logic_name => 'gzip_non_chromosomal',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd => q(gzip -c #ftp_dumps_dir#/non_chromosomal/*gff >#output_dir#/#species#.#assembly#.nonchromosomal.motif_features.gff.gz  )
            },
            -rc_name    => '32Gb_job'
        },
        {
            -logic_name => 'tabix',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd => q(tabix #output_dir#/#species#.#assembly#.motif_features.gff.gz)
            },
            -rc_name    => '32Gb_job'
        }

    ];
}

1;