=pod 

=head1 NAME

    Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::ChIPSeqCleanup

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 LICENSE

    Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
    Copyright [2016-2017] EMBL-European Bioinformatics Institute

    Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

         http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software distributed under the License
    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and limitations under the License.

=head1 CONTACT

=cut

package Bio::EnsEMBL::Funcgen::PipeConfig::RegulatoryBuild::MapStableIds;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base 'Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::Base';

sub pipeline_analyses {
    my ($self) = @_;
    return [
        {   -logic_name => 'start_regulatory_build_stable_id_mapping',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into => { 
              MAIN => 'regulatory_build_stable_id_mapping_job_factory',
            },
        },
        {   -logic_name => 'regulatory_build_stable_id_mapping_job_factory',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -parameters => {
                db_conn    => 'funcgen:#species#',
            },
            -flow_into => {
                MAIN => 'export_regulatory_features_to_bed',
            },
        },
        {   -logic_name => 'export_regulatory_features_to_bed',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
              cmd => qq( 
                export_regulatory_features_to_bed.pl \
                  -url              '#db_url_funcgen#' \
                  -species          #species# \
                  -stable_id_prefix #stable_id_prefix# \
                  -outfile          #regulatory_features_bed_file#
              )
            },
            -flow_into => {
                MAIN => 'export_regulatory_features_to_bed_old',
            },
        },
        {   -logic_name => 'export_regulatory_features_to_bed_old',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
              cmd => qq( 
                export_regulatory_features_to_bed.pl \
                  -url              '#db_url_funcgen_old#' \
                  -species          #species# \
                  -stable_id_prefix #stable_id_prefix# \
                  -outfile          #regulatory_features_previous_version_bed_file#
              )
            },
            -flow_into => {
                MAIN => 'sort_bed_new',
            },
        },
        {   -logic_name => 'sort_bed_new',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
              cmd => qq( 
                bedSort #regulatory_features_bed_file# #regulatory_features_bed_file#
              )
            },
            -flow_into => {
                MAIN => 'sort_bed_old',
            },
        },
        {   -logic_name => 'sort_bed_old',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
              cmd => qq( 
                bedSort #regulatory_features_previous_version_bed_file# #regulatory_features_previous_version_bed_file#
              )
            },
            -flow_into => {
                MAIN => 'compute_overlaps',
            },
        },
        {   -logic_name => 'compute_overlaps',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
              cmd => qq( 
                bedtools intersect -a #regulatory_features_previous_version_bed_file# -b #regulatory_features_bed_file# -wo > #overlaps_bed_file#
              )
            },
            -flow_into => {
                MAIN => 'map_stable_ids',
            },
        },
        {   -logic_name => 'map_stable_ids',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
              cmd => qq( 
                generate_stable_ids_from_overlaps_using_length.pl \
                    --all_overlaps               #overlaps_bed_file# \
                    --source_regulatory_features #regulatory_features_previous_version_bed_file# \
                    --target_regulatory_features #regulatory_features_bed_file# \
                    --stable_id_prefix           #stable_id_prefix# \
                    --outfile                    #stable_id_mapping_file#
              )
            },
            -flow_into => {
                MAIN => 'remove_old_stable_ids_in_db',
            },
        },
        {   -logic_name => 'remove_old_stable_ids_in_db',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
              cmd => qq( 
                set_stable_ids_to_null_in_db.pl \
                    --url     '#db_url_funcgen#' \
                    --species #species#
              )
            },
            -flow_into => {
                MAIN => 'load_stable_ids',
            },
        },
        {   -logic_name => 'load_stable_ids',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
              cmd => qq( 
                update_stable_ids_in_db.pl \
                    --url          '#db_url_funcgen#' \
                    --species      #species# \
                    --mapping_file #stable_id_mapping_file#
              )
            },
            -flow_into => {
                MAIN => 'regulatory_build_stable_id_mapping_done',
            },
        },
        {   -logic_name => 'regulatory_build_stable_id_mapping_done',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
    ];
}

1;
