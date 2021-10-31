=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::PipeConfig::MotifMapping::MOODS_conf

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

package Bio::EnsEMBL::Funcgen::PipeConfig::MotifMapping::MOODS_conf;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base ('Bio::EnsEMBL::Funcgen::PipeConfig::MotifMapping::Base_conf');

sub pipeline_analyses {
    my $self = shift;

    return [
        {
            -logic_name => 'start_MOODS',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into  => {
                '1->A' => [ 'dump_binding_matrix_frequencies_to_disk' ],
                'A->1' => 'matrix_factory_1'
            },

        },
        {
            -logic_name => 'dump_binding_matrix_frequencies_to_disk',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::'
                            . 'DumpBindingMatrixFrequenciesToDisk',
        },
        {
            -logic_name => 'matrix_factory_1',
            # -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::MotifMapping::BindingMatrixJobFactory',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -parameters => { inputquery   => q(SELECT stable_id FROM binding_matrix),
                             db_conn      => 'funcgen:#species#',
                             column_names => [ 'matrix' ]
            },
            -flow_into  => {
                2 => { 'run_MOODS' => INPUT_PLUS() },
            }
        },
        {
            -logic_name => 'run_MOODS',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => { cmd =>
                             'python2.7 /nfs/production/flicek/ensembl/regulation/pipelines/motif_features/MOODS-python-1.9.4.1/scripts/moods-dna.py '
                                 . '-m #matrices_dir#/#matrix#.pfm '
                                 . '-s #genome# '
                                 . '-p 0.01 '
                                 . '>#MOODS_raw_dir#/#matrix#.out'
            },
            -rc_name    => '8Gb_job',
            -flow_into  => {
                1 => { 'convert_2_sorted_BED' => INPUT_PLUS() },
            }
        },
        {
            -logic_name => 'convert_2_sorted_BED',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd =>
                q(awk 'BEGIN{FS=","; OFS="\t"}; { split($1, contig_info, " "); len=length($6); print contig_info[1],$3+1,$3+len,$5"_"contig_info[3],"1000",$4}' )
                    . ' #MOODS_raw_dir#/#matrix#.out >#MOODS_bed_dir#/#matrix#.bed '
                    . ' && sort -k1,1 -k2,2n -k3,3n #MOODS_bed_dir#/#matrix#.bed >#MOODS_bed_dir#/#matrix#.sorted.bed '
                    . ' && rm #MOODS_bed_dir#/#matrix#.bed'
            },
            -rc_name    => '16Gb_job',
            -max_retry_count => 1,
            -flow_into       => {
                MEMLIMIT => 'convert_2_sorted_BED_highmem'
            }
        },
        {
            -logic_name => 'convert_2_sorted_BED_highmem',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd =>
                q(awk 'BEGIN{FS=","; OFS="\t"}; { split($1, contig_info, " "); len=length($6); print contig_info[1],$3+1,$3+len,$5"_"contig_info[3],"1000",$4}' )
                    . ' #MOODS_raw_dir#/#matrix#.out >#MOODS_bed_dir#/#matrix#.bed '
                    . ' && sort -k1,1 -k2,2n -k3,3n #MOODS_bed_dir#/#matrix#.bed >#MOODS_bed_dir#/#matrix#.sorted.bed '
                    . ' && rm #MOODS_bed_dir#/#matrix#.bed'
            },
            -rc_name    => '32Gb_job',
            -max_retry_count => 1,
        }

    ];
}

1;
