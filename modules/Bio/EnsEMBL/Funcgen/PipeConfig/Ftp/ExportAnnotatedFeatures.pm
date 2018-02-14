package Bio::EnsEMBL::Funcgen::PipeConfig::Ftp::ExportAnnotatedFeatures;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Funcgen::PipeConfig::Ftp::Base';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use Hash::Util qw( lock_hash );

sub pipeline_analyses {
    my $self = shift;
    
    my $data_freeze_date = $self->o('data_freeze_date');
    
    my $ftp_layout_configuration = {
      peaks_gff_file_dir    => '#ftp_base_dir#/#species#/Peaks/#epigenome_production_name#/#feature_type_name#',
      peaks_bed_file_dir    => '#ftp_base_dir#/#species#/Peaks/#epigenome_production_name#/#feature_type_name#',
      peaks_bigbed_file_dir => '#ftp_base_dir#/#species#/Peaks/#epigenome_production_name#/#feature_type_name#',
      
      peaks_gff_file_base_name    => "#species#.#assembly#.#epigenome_production_name#.#feature_type_name#.#analysis_logic_name#.peaks.${data_freeze_date}.gff.gz",
      peaks_bed_file_base_name    => "#species#.#assembly#.#epigenome_production_name#.#feature_type_name#.#analysis_logic_name#.peaks.${data_freeze_date}.bed.gz",
      peaks_bigbed_file_base_name => "#species#.#assembly#.#epigenome_production_name#.#feature_type_name#.#analysis_logic_name#.peaks.${data_freeze_date}.bb",
    };

    lock_hash(%$ftp_layout_configuration);

    return [
        {   -logic_name  => 'start_export',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => 'mk_temp_dir',
        },
        {   -logic_name  => 'mk_temp_dir',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
              cmd => 'mkdir -p #tempdir_ftp#',
            },
            -flow_into   => {
               MAIN => 'create_chromosome_length_file_for_bigbed_conversion',
            },
        },
        {   -logic_name  => 'create_chromosome_length_file_for_bigbed_conversion',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
              cmd => 'create_chromosome_length_file_for_bigbed_conversion.pl --registry #reg_conf# --species #species# > #tempdir_ftp#/#species#.ucsc.sizes',
            },
            -flow_into   => {
               MAIN => 'for_each_epigenome_and_feature_type_having_peaks',
            },
        },
        {   -logic_name  => 'for_each_epigenome_and_feature_type_having_peaks',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::Ftp::JobFactoryAnnotatedFeaturesUsingIdLists',
            -flow_into   => {
               '2->A' => 'make_temp_dir',
               'A->2' => 'rm_peaks_temp_dir',
            },
        },
        {   -logic_name  => 'make_temp_dir',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
              cmd => 'mkdir -p #tempdir_ftp#/#epigenome_production_name#/#feature_type_name#',
            },
            -flow_into   => {
               MAIN => 'print_annotated_feature_ids',
            },
        },
        {   -logic_name  => 'print_annotated_feature_ids',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -rc_name           => '2Gb_job',
            -parameters  => {
              cmd => '
                print_annotated_feature_ids.pl  \
                  --registry #reg_conf#  \
                  --species #species#  \
                  --epigenome_production_name #epigenome_production_name#  \
                  --feature_type_name #feature_type_name#  \
                  > #tempdir_ftp#/#epigenome_production_name#/#feature_type_name#/#epigenome_production_name#_#feature_type_name#_ids.txt
              ',
            },
            -flow_into   => {
               MAIN => 'export_annotated_features',
            },
        },
        {   -logic_name  => 'export_annotated_features',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -analysis_capacity => 10,
            -batch_size        => 1,
            -rc_name           => '2Gb_job',
            -parameters  => {
              cmd => '
                export_annotated_features.pl \
                  --gff_file #tempdir_ftp#/#epigenome_production_name#/#feature_type_name#/#epigenome_production_name#_#feature_type_name#.peaks.unsorted.gff  \
                  --bed_file #tempdir_ftp#/#epigenome_production_name#/#feature_type_name#/#epigenome_production_name#_#feature_type_name#.peaks.unsorted.bed  \
                  --registry #reg_conf#  \
                  --species #species#  \
                  --ids #tempdir_ftp#/#epigenome_production_name#/#feature_type_name#/#epigenome_production_name#_#feature_type_name#_ids.txt
              ',
            },
            -flow_into   => {
              MAIN => [
                'sort_peaks_gff',
                'sort_peaks_bed',
              ],
            },
        },
        {   -logic_name  => 'sort_peaks_gff',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'sort -k1,1 -k4,4n -o #tempdir_ftp#/#epigenome_production_name#/#feature_type_name#/#epigenome_production_name#_#feature_type_name#.peaks.gff #tempdir_ftp#/#epigenome_production_name#/#feature_type_name#/#epigenome_production_name#_#feature_type_name#.peaks.unsorted.gff ; sleep 10',
            },
            -flow_into   => {
               MAIN => 'gzip_peaks_gff',
            },
            -rc_name    => '2Gb_job',
        },
        {   -logic_name  => 'gzip_peaks_gff',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'gzip #tempdir_ftp#/#epigenome_production_name#/#feature_type_name#/#epigenome_production_name#_#feature_type_name#.peaks.gff',
            },
            -flow_into   => {
               MAIN => 'mk_peaks_gff_ftp_dir',
            },
        },
        {   -logic_name  => 'mk_peaks_gff_ftp_dir',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'mkdir -p ' . $ftp_layout_configuration->{peaks_gff_file_dir},
            },
            -flow_into   => {
               MAIN => 'mv_peaks_gff_to_ftp_dir',
            },
        },
        {   -logic_name  => 'mv_peaks_gff_to_ftp_dir',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'mv #tempdir_ftp#/#epigenome_production_name#/#feature_type_name#/#epigenome_production_name#_#feature_type_name#.peaks.gff.gz ' . $ftp_layout_configuration->{peaks_gff_file_dir} . '/' . $ftp_layout_configuration->{peaks_gff_file_base_name},
            },
        },

        {   -logic_name  => 'sort_peaks_bed',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => '
                  bedSort  \
                    #tempdir_ftp#/#epigenome_production_name#/#feature_type_name#/#epigenome_production_name#_#feature_type_name#.peaks.unsorted.bed  \
                    #tempdir_ftp#/#epigenome_production_name#/#feature_type_name#/#epigenome_production_name#_#feature_type_name#.peaks.bed ; sleep 10
                ',
            },
            -flow_into   => {
               MAIN => 'convert_to_bigbed',
            },
            -rc_name    => '2Gb_job',
        },
        {   -logic_name  => 'convert_to_bigbed',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => '
                  bedToBigBed  \
                    #tempdir_ftp#/#epigenome_production_name#/#feature_type_name#/#epigenome_production_name#_#feature_type_name#.peaks.bed  \
                    #tempdir_ftp#/#species#.ucsc.sizes  \
                    #tempdir_ftp#/#epigenome_production_name#/#feature_type_name#/#epigenome_production_name#_#feature_type_name#.peaks.bb
                ',
            },
            -flow_into   => {
               MAIN => 'gzip_peaks_bed',
            },
            -rc_name    => '2Gb_job',
        },
        {   -logic_name  => 'gzip_peaks_bed',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'gzip #tempdir_ftp#/#epigenome_production_name#/#feature_type_name#/#epigenome_production_name#_#feature_type_name#.peaks.bed',
            },
            -flow_into   => {
               MAIN => [
                'mk_peaks_bed_to_ftp_dir',
                'mk_peaks_bigbed_to_ftp_dir',
               ]
            },
        },
        {   -logic_name  => 'mk_peaks_bed_to_ftp_dir',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'mkdir -p ' . $ftp_layout_configuration->{peaks_bed_file_dir},
            },
            -flow_into   => {
               MAIN => 'mv_peaks_bed_to_ftp',
            },
        },
        {   -logic_name  => 'mk_peaks_bigbed_to_ftp_dir',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'mkdir -p ' . $ftp_layout_configuration->{peaks_bigbed_file_dir},
            },
            -flow_into   => {
               MAIN => 'mv_annotated_feature_bigbed_to_ftp'
            },
        },
        {   -logic_name  => 'mv_peaks_bed_to_ftp',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'mv #tempdir_ftp#/#epigenome_production_name#/#feature_type_name#/#epigenome_production_name#_#feature_type_name#.peaks.bed.gz ' . $ftp_layout_configuration->{peaks_bed_file_dir} . '/' . $ftp_layout_configuration->{peaks_bed_file_base_name},
            },
        },
        {   -logic_name  => 'mv_annotated_feature_bigbed_to_ftp',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'mv #tempdir_ftp#/#epigenome_production_name#/#feature_type_name#/#epigenome_production_name#_#feature_type_name#.peaks.bb ' . $ftp_layout_configuration->{peaks_bigbed_file_dir} . '/' . $ftp_layout_configuration->{peaks_bigbed_file_base_name},
            },
        },
        {
          -logic_name => 'rm_peaks_temp_dir',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => 'rm -rf #tempdir_ftp#/#epigenome_production_name#/#feature_type_name#',
          },
        },
    ]
}

1;
