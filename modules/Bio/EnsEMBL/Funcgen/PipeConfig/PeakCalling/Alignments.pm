package Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::Alignments;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::Base';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub pipeline_analyses {
    my $self = shift;

    return [
        {   -logic_name  => 'start_alignments',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               MAIN => 'truncate_alignment_tables',
            },
        },
        {
            -logic_name  => 'truncate_alignment_tables',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
            -parameters => {
                sql     => [
                    "truncate peak;",
                    "truncate peak_calling;",
                    "truncate peak_calling_statistic",
                    "truncate alignment;",
                    "delete from data_file where table_name = 'alignment';",
                    "truncate alignment_qc_flagstats;",
                    "truncate alignment_read_file;",
                    "truncate chance;",
                    "truncate fastqc;",
                    "truncate frip;",
                    "truncate idr;",
                    "truncate phantom_peak;",
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
                    "truncate segmentation_cell_tables;",
                ],
                db_conn => 'funcgen:#species#',
            },
            -flow_into   => {
               MAIN => 'seed_control_experiments'
            },
        },

        {   -logic_name  => 'seed_control_experiments',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::FindControlExperiments',
            -flow_into   => {
               'A->2' => 'start_align_signals',
               '3->A' => 'start_align_controls',
               4      => 'start_align_signals',
            },
        },
        {
            -logic_name  => 'start_align_controls',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        {
            -logic_name  => 'start_align_signals',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
    ]
}

1;
