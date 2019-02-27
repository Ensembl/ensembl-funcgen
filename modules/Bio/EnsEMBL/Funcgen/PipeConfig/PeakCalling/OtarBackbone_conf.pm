package Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::OtarBackbone_conf;

use strict;
use warnings;
use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub beekeeper_extra_cmdline_options {
  my $self = shift;
  return '-reg_conf ' . $self->o('reg_conf') . ' -keep_alive -can_respecialize 1 -sleep 0.2';
}

sub default_options {
  my $self = shift;
  
  return {
      %{$self->SUPER::default_options},
      pipeline_name => 'otar',
   };
}

sub pipeline_wide_parameters {
  my $self = shift;
  
  return {
    %{$self->SUPER::pipeline_wide_parameters},
    pipeline_name => $self->o('pipeline_name'),
  };
}

sub pipeline_analyses {
    my $self = shift;

    return [
        {   -logic_name  => 'start',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               MAIN => 'backbone_fire_pre_pipeline_checks',
            },
        },
        {   -logic_name  => 'backbone_fire_pre_pipeline_checks',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_pre_pipeline_checks',
               'A->1' => 'backbone_fire_populate_read_file_stats'
            },
        },
        {
            -logic_name => 'start_pre_pipeline_checks',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },

        {   -logic_name  => 'backbone_fire_populate_read_file_stats',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_populate_read_file_stats',
               'A->1' => 'backbone_fire_fastqc'
            },
        },
        {
            -logic_name => 'start_populate_read_file_stats',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },

        {   -logic_name  => 'backbone_fire_fastqc',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_fastqc',
               'A->1' => 'backbone_fire_fastqc_report'
            },
        },
        {
            -logic_name => 'start_fastqc',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },



        {   -logic_name  => 'backbone_fire_fastqc_report',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_fastqc_report',
               'A->1' => 'backbone_fire_alignments'
            },
        },
        {
            -logic_name => 'start_fastqc_report',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },







        {   -logic_name  => 'backbone_fire_alignments',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_alignments',
               'A->1' => 'backbone_fire_alignment_hc'
            },
        },
        {
            -logic_name => 'start_alignments',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },

        {   -logic_name  => 'backbone_fire_alignment_hc',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_alignment_hc',
               'A->1' => 'backbone_fire_write_bigwig'
            },
        },
        {
            -logic_name => 'start_alignment_hc',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        {   -logic_name  => 'backbone_fire_write_bigwig',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_write_bigwig',
               'A->1' => 'backbone_fire_alignment_qc'
            },
        },
        {
            -logic_name => 'start_write_bigwig',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },

        {   -logic_name  => 'backbone_fire_alignment_qc',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_alignment_qc',
               'A->1' => 'backbone_fire_peak_calling'
            },
        },
        {
            -logic_name => 'start_alignment_qc',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },

        {   -logic_name  => 'backbone_fire_peak_calling',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_peak_calling',
               'A->1' => 'backbone_fire_peak_calling_hc'
            },
        },
        {
            -logic_name => 'start_peak_calling',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        {   -logic_name  => 'backbone_fire_peak_calling_hc',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_peak_calling_hc',
               'A->1' => 'backbone_fire_frip'
            },
        },
        {
            -logic_name  => 'start_peak_calling_hc',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        {   -logic_name  => 'backbone_fire_frip',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_frip',
               'A->1' => 'backbone_fire_cleanup'
            },
        },
        {
            -logic_name  => 'start_frip',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        {   -logic_name  => 'backbone_fire_cleanup',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_cleanup',
               'A->1' => 'backbone_fire_quality_check_reports'
            },
        },
        {
            -logic_name  => 'start_cleanup',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        {   -logic_name  => 'backbone_fire_quality_check_reports',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_quality_check_reports',
               'A->1' => 'backbone_fire_segmentation'
            },
        },
        {
            -logic_name  => 'start_quality_check_reports',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        {   -logic_name  => 'backbone_fire_segmentation',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_segmentation',
               'A->1' => 'backbone_fire_segmentation_statistics'
            },
        },
        {
            -logic_name  => 'start_segmentation',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },

        {   -logic_name  => 'backbone_fire_segmentation_statistics',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_segmentation_statistics',
               'A->1' => 'backbone_fire_regulatory_build_hc'
            },
        },
        {
            -logic_name  => 'start_segmentation_statistics',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },

        {   -logic_name  => 'backbone_fire_regulatory_build_hc',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_regulatory_build_hc',
               'A->1' => 'backbone_fire_regulatory_build_statistics'
            },
        },
        {
            -logic_name  => 'start_regulatory_build_hc',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        
        {   -logic_name  => 'backbone_fire_regulatory_build_statistics',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_regulatory_build_statistics',
               'A->1' => 'backbone_fire_regulatory_build_stable_id_mapping'
            },
        },
        {
            -logic_name  => 'start_regulatory_build_statistics',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        
        {   -logic_name  => 'backbone_fire_regulatory_build_stable_id_mapping',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_regulatory_build_stable_id_mapping',
               'A->1' => 'backbone_fire_stable_id_mapping_hc'
            },
        },
        {
            -logic_name  => 'start_regulatory_build_stable_id_mapping',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        {   -logic_name  => 'backbone_fire_stable_id_mapping_hc',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_stable_id_mapping_hc',
               'A->1' => 'backbone_fire_ftp_export'
            },
        },
        {
            -logic_name  => 'start_stable_id_mapping_hc',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
        {   -logic_name  => 'backbone_fire_ftp_export',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_ftp_export',
               'A->1' => 'backbone_pipeline_finished'
            },
        },
        {
            -logic_name  => 'start_ftp_export',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into  => {
              MAIN => 'backbone_fire_exports',
             }
        },
        {
            -logic_name  => 'backbone_fire_exports',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },

        {   -logic_name => 'backbone_pipeline_finished',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        }
    ]
}

1;
