=pod 

=head1 NAME

    Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::ChIPSeqCleanup

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 LICENSE

    Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
    Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

sub pipeline_wide_parameters {
  my $self = shift;
  
  return {
    %{$self->SUPER::pipeline_wide_parameters},
    pipeline_name            => $self->o('pipeline_name'),
    tempdir                  => $self->o('tempdir'),
    tempdir_regulatory_build => $self->o('tempdir') . '/regulatory_build',
    data_root_dir            => $self->o('data_root_dir'),
    reg_conf                 => $self->o('reg_conf'),
  };
}

sub pipeline_analyses {
    my ($self) = @_;
    
    my $previous_version_species_suffix = '_previous_version';
    
    return [
        {   -logic_name => 'start_regulatory_build_stable_id_mapping',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into => { 
              MAIN => 'regulatory_build_stable_id_mapping_job_factory',
            },
        },
        {   -logic_name => 'regulatory_build_stable_id_mapping_job_factory',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::RegulatoryBuild::StableIdMappingJobFactory',
            -parameters => {
                db_conn => 'funcgen:#species#',
                tempdir => '#tempdir_regulatory_build#/#species#/stable_id_mapping',
            },
            -flow_into => {
                '2->A' => [
                    'export_regulatory_features_to_bed',
                    'export_regulatory_features_to_bed_old',
                  ],
                'A->2' => [ 'compute_overlaps' ]
            },
        },
        {   -logic_name => 'export_regulatory_features_to_bed',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
              cmd => q( 
                export_regulatory_features_to_bed.pl \
                  -registry         #reg_conf# \
                  -species          #species# \
                  -stable_id_prefix #stable_id_prefix# \
                  -outfile          #regulatory_features_bed_file#
              )
            },
            -flow_into => {
                MAIN => 'sort_bed_new',
            },
        },
        {   -logic_name => 'export_regulatory_features_to_bed_old',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
              cmd => qq( 
                export_regulatory_features_to_bed.pl \\
                  -registry         #reg_conf# \\
                  -species          #species#$previous_version_species_suffix \\
                  -stable_id_prefix #stable_id_prefix# \\
                  -outfile          #regulatory_features_previous_version_bed_file#
              )
            },
            -flow_into => {
                MAIN => 'sort_bed_old',
            },
        },
        {   -logic_name => 'sort_bed_new',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
              cmd => qq( 
                bedSort #regulatory_features_bed_file# #regulatory_features_bed_file#
              )
            },
        },
        {   -logic_name => 'sort_bed_old',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
              cmd => qq( 
                bedSort #regulatory_features_previous_version_bed_file# #regulatory_features_previous_version_bed_file#
              )
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
              cmd => q( 
                generate_stable_ids_from_overlaps_using_length.pl \
                    --all_overlaps               #overlaps_bed_file# \
                    --source_regulatory_features #regulatory_features_previous_version_bed_file# \
                    --target_regulatory_features #regulatory_features_bed_file# \
                    --stable_id_prefix           #stable_id_prefix# \
                    --outfile                    #stable_id_mapping_file# \
                    --mapping_report             #mapping_report#
              )
            },
            -rc_name    => '32Gb_job',
            -flow_into => {
                MAIN => [
                    'remove_old_stable_ids_in_db',
                    'store_stable_id_mapping_statistics',
                ],
            },
        },
        {   -logic_name => 'store_stable_id_mapping_statistics',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
              cmd => q( 
                store_stable_id_mapping_statistics.pl \
                    --species             #species# \
                    --registry            #reg_conf# \
                    --mapping_report_file #mapping_report#
              )
            },
            -flow_into => {
                MAIN => 'generate_regulatory_build_report_with_stable_id_mapping',
            },
        },
      {   -logic_name => 'generate_regulatory_build_report_with_stable_id_mapping',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 
                q(
                  generate_regulatory_build_report.pl  \
                      --species     #species# \
                      --registry    #reg_conf# \
                      --output_file #reports_dir#/#species#/regulatory_build.html
                )
          },
        },
        {   -logic_name => 'remove_old_stable_ids_in_db',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
              cmd => qq( 
                set_stable_ids_to_null_in_db.pl -registry #reg_conf# --species #species#
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
                update_stable_ids_in_db.pl -registry #reg_conf# --species #species# --mapping_file #stable_id_mapping_file#
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
