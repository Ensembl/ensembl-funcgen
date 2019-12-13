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

package Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::PeakCallingHc;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base 'Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::Base';

sub pipeline_analyses {
    my ($self) = @_;
    return [
        {   -logic_name => 'start_peak_calling_hc',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into => { 
              MAIN => 'truncate_peak_calling_statistics',
            },
        },
        {   -logic_name => 'truncate_peak_calling_statistics',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
            -parameters => {
                db_conn => 'funcgen:#species#',
                sql     => [
                    q~
                        truncate peak_calling_statistic;
                    ~,
                ]
            },
            -flow_into => { 
              MAIN => 'compute_peak_calling_statistics',
            },
        },
        {
            -logic_name  => 'compute_peak_calling_statistics',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
            -parameters => {
                sql     => [
                    q~
                        insert into peak_calling_statistic (peak_calling_id, epigenome_id, feature_type_id, statistic, value)
                        select
                            peak_calling_id,
                            null,
                            feature_type_id,
                            'total_length',
                            sum(peak.seq_region_end - peak.seq_region_start + 1) as total_length
                        from
                            peak_calling
                            left join peak using (peak_calling_id)
                        group by
                            peak_calling_id
                        ;
                    ~,
                    q~
                        insert into peak_calling_statistic (peak_calling_id, epigenome_id, feature_type_id, statistic, value)
                        select
                            peak_calling_id,
                            null,
                            feature_type_id,
                            'num_peaks',
                            count(peak.peak_id) as num_peaks
                        from
                            peak_calling
                            left join peak using (peak_calling_id)
                        group by
                            peak_calling_id
                        ;
                    ~,
                    q~
                        insert into peak_calling_statistic (peak_calling_id, epigenome_id, feature_type_id, statistic, value)
                        select
                            peak_calling_id,
                            null,
                            feature_type_id,
                            'average_length',
                            avg(peak.seq_region_end - peak.seq_region_start + 1) as average_length
                        from
                            peak_calling
                            left join peak using (peak_calling_id)
                        group by
                            peak_calling_id
                        ;
                    ~,
                ],
                db_conn => 'funcgen:#species#',
            },
            -flow_into => { 
              MAIN => 'peak_calling_statistics_job_factory',
            },
        },
        {   -logic_name => 'peak_calling_statistics_job_factory',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -parameters => { 
                inputlist => [
                    qw(
                        CTCF
                        DNase1
                        H3K27ac
                        H3K27me3
                        H3K36me3
                        H3K4me1
                        H3K4me2
                        H3K4me3
                        H3K9ac
                        H3K9me3
                    )
                ],
            },
            -flow_into => {
              '2->A'    => { 'compute_overall_peak_calling_statistics' => INPUT_PLUS },
              'A->MAIN' => 'generate_peak_calling_report',
            },
        },
        {   -logic_name => 'compute_overall_peak_calling_statistics',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -rc_name    => '4Gb_job',
            -parameters => { 
              cmd => 
                  qq( compute_peak_calling_statistics.pl                           )
                . qq(   --species      #species#                                   )
                . qq(   --registry     #reg_conf#                                  )
                . qq(   --feature_type #_0#                                        )
                . qq(   --tempdir      #tempdir_peak_calling#/#species#/peak_stats )
            },
        },
        {   -logic_name => 'generate_peak_calling_report',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => { 
              cmd => qq( generate_peak_calling_report.pl )
                . qq( --species          #species#         )
                . qq( --registry         #reg_conf#        )
                . qq( --output_directory #reports_dir#/#species# )
            },
            -flow_into => {
              MAIN     => 'hc_only_bigwig_files_for_complete_and_deduplicated_alignments',
            },
        },
        {
            -logic_name  => 'hc_only_bigwig_files_for_complete_and_deduplicated_alignments',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SqlHealthcheck',
            -parameters => {
              db_conn       => 'funcgen:#species#',
              description   => 'Bigwig files should only exist for complete alignments.',
              query         => "select * from alignment where bigwig_file_id is null and is_complete = true and has_duplicates = false",
              expected_size => '0'
            },
          -flow_into => {
              MAIN => 'hc_peak_calls_for_all_signal_experiments_available',
          },
        },
        {
            -logic_name  => 'hc_peak_calls_for_all_signal_experiments_available',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SqlHealthcheck',
            -parameters => {
              db_conn       => 'funcgen:#species#',
              description   => 'Bigwig files should only exist for complete alignments.',
              query         => "select experiment.experiment_id from experiment left join peak_calling using (experiment_id) where peak_calling.experiment_id is null and experiment.is_control = 0",
              expected_size => '0'
            },
            -failed_job_tolerance => 100,
          -flow_into => {
              MAIN => 'peak_calling_hc_done',
          },
        },
        {   -logic_name => 'peak_calling_hc_done',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
    ];
}

1;
