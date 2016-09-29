package Bio::EnsEMBL::Funcgen::HiveConfig::Ftp::ExportRegulatoryFeatures;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Funcgen::HiveConfig::Ftp::Base';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub pipeline_analyses {
    my $self = shift;
    
    my $data_freeze_date = $self->o('data_freeze_date');

    my $ftp_layout_configuration = {
      regulatory_activities_gff_file_dir       => '#ftp_base_dir#/#species#/RegulatoryFeatureActivity/#epigenome_production_name#',
      regulatory_activities_gff_file_base_name => "#species#.#assembly#.#epigenome_production_name#.#analysis_logic_name#.regulatory_activity.${data_freeze_date}.gff",

      regulatory_features_gff_file_dir       => '#ftp_base_dir#/#species#',
      regulatory_features_gff_file_base_name => "#species#.#assembly#.Regulatory_Build.regulatory_features.${data_freeze_date}.gff",
    };

    return [
        {   -logic_name  => 'start_export',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               MAIN => 'job_factory_regulatory',
            },
        },
        {   -logic_name  => 'job_factory_regulatory',
            -module      => 'Bio::EnsEMBL::Funcgen::Hive::Ftp::JobFactoryQualityChecks',
            -flow_into   => {
               MAIN => { 
                  'dbconn_for_species'         => INPUT_PLUS(),
                  'export_regulatory_features' => INPUT_PLUS(),
                }
             }
        },
        {   -logic_name  => 'dbconn_for_species',
            -module      => 'Bio::EnsEMBL::Funcgen::Hive::Ftp::DbconnForSpecies',
            -parameters  => {
                group => 'funcgen',
            },
            -flow_into   => {
               2 => {
                'job_factory_regulatory_activities' => INPUT_PLUS(),
               }
            },
        },
        {   -logic_name  => 'job_factory_regulatory_activities',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -parameters  => { 
	      db_conn    => '#url#',
	      inputquery => 'select epigenome.name as epigenome_name, epigenome.production_name as epigenome_production_name, analysis.logic_name as analysis_logic_name from regulatory_build join regulatory_build_epigenome using (regulatory_build_id) join epigenome using (epigenome_id) join analysis using (analysis_id) where regulatory_build.is_current=1',
            },
            -flow_into   => {
               2 => { 'export_regulatory_activities', INPUT_PLUS() },
            },
        },
        {   -logic_name  => 'export_regulatory_activities',
	    -analysis_capacity => 20,
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'export_regulatory_features.pl --epigenome_name "#epigenome_name#" --output_file '
                  . $ftp_layout_configuration->{regulatory_activities_gff_file_dir} . '/' . $ftp_layout_configuration->{regulatory_activities_gff_file_base_name} 
                  . ' --registry #reg_conf# --species #species#',
            },
            -flow_into   => {
               MAIN => 'gzip_regulatory_activities'
            },
        },
        {   -logic_name  => 'gzip_regulatory_activities',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'gzip ' . $ftp_layout_configuration->{regulatory_activities_gff_file_dir} . '/' . $ftp_layout_configuration->{regulatory_activities_gff_file_base_name},
            },
        },
        {   -logic_name  => 'export_regulatory_features',
	    -analysis_capacity => 20,
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'export_regulatory_features.pl --only_summary --output_file '
                  . $ftp_layout_configuration->{regulatory_features_gff_file_dir} . '/' . $ftp_layout_configuration->{regulatory_features_gff_file_base_name} 
                  . ' --registry #reg_conf# --species #species#',
            },
            -flow_into   => {
               MAIN => 'gzip_regulatory_features'
            },
        },
        {   -logic_name  => 'gzip_regulatory_features',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'gzip ' . $ftp_layout_configuration->{regulatory_features_gff_file_dir} . '/' . $ftp_layout_configuration->{regulatory_features_gff_file_base_name},
            },
        },

    ]
}

1;
