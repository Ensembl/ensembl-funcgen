package Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::PeakCalling;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::Base';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub pipeline_analyses {
    my $self = shift;

    return [
        {
            -logic_name  => 'start_peak_calling',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               'MAIN->A' => 'truncate_peak_calling_tables',
               'A->MAIN' => 'presort_peak_table',
            },
        },
        {
            -logic_name => 'presort_peak_table',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd => qq( presort_peak_table.pl    )
                     . qq(   --species  #species#   )
                     . qq(   --registry #reg_conf#  )
            },
            -flow_into   => {
               MAIN => 'meta_coord_for_peaks',
            },
        },
        {
            -logic_name => 'meta_coord_for_peaks',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd => qq( populate_meta_coord.pl    )
                     . qq(   --species  #species#    )
                     . qq(   --registry #reg_conf#   )
            },
        },
        {
            -logic_name  => 'truncate_peak_calling_tables',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
            -parameters => {
                sql     => [
                    "truncate peak;",
                    "truncate peak_calling;",
                    "truncate peak_calling_statistic",
                    "truncate alignment_qc_flagstats;",
                    "truncate frip;",
                    "truncate idr;",
                    "truncate regulatory_build;",
                    "truncate regulatory_feature;",
                    "truncate regulatory_activity;",
                    "truncate regulatory_build_epigenome;",
                    "truncate regulatory_evidence;",
                    "truncate regulatory_build_statistic;",
                    "truncate segmentation_file;",
                    "delete from data_file where table_name = 'segmentation_file';",
                    "truncate segmentation_state_assignment;",
                    "truncate segmentation_state_emission;",
                ],
                db_conn => 'funcgen:#species#',
            },
            -flow_into   => {
               MAIN => 'seed_peak_calling_jobs_from_list',
            },
        },
        {
            -logic_name  => 'seed_peak_calling_jobs_from_list',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::SeedAllExecutionPlans',
            -flow_into   => {
               2 => 'backbone_fire_convert_signal_to_bed',
            },
        },

        {   -logic_name  => 'backbone_fire_convert_signal_to_bed',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_convert_to_bed',
               'A->1' => 'backbone_fire_idr'
            },
        },
        {
            -logic_name  => 'start_convert_to_bed',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },

        {   -logic_name  => 'backbone_fire_idr',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_idr',
               'A->1' => 'backbone_fire_call_peaks'
            },
        },
        {
            -logic_name  => 'start_idr',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },

        {   -logic_name  => 'backbone_fire_call_peaks',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_call_peaks',
               'A->1' => 'backbone_chipseq_finished'
            },
        },
        {
            -logic_name  => 'start_call_peaks',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },

        {   -logic_name => 'backbone_chipseq_finished',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        }
    ]
}

1;
