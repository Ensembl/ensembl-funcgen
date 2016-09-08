package Bio::EnsEMBL::Funcgen::HiveConfig::Ftp::ExportRegulatoryFeatures;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Hive::PipeConfig::EnsemblGeneric_conf';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
      %{$self->SUPER::pipeline_wide_parameters},

      ftp_base_dir  => $self->o('ftp_base_dir'),
      reg_conf      => $self->o('reg_conf'),
    };
}

sub beekeeper_extra_cmdline_options {
    my ($self) = @_;
    return '-reg_conf ' . $self->o('reg_conf') . ' -keep_alive -can_respecialize 1';
}

sub pipeline_analyses {
    my $self = shift;

    return [
        {   -logic_name  => 'start_export',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               'MAIN' => [ 'dbconn_for_species', 'export_regulatory_features' ]
            },
        },
        {   -logic_name  => 'dbconn_for_species',
            -module      => 'Bio::EnsEMBL::Funcgen::Hive::Ftp::DbconnForSpecies',
            -parameters  => { 
	      group => 'funcgen'
            },
            -flow_into   => {
               2 => { 'job_factory_regulatory_activities', INPUT_PLUS() },
            },
        },
        {   -logic_name  => 'job_factory_regulatory_activities',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -parameters  => { 
	      db_conn    => '#url#',
	      inputquery => 'select epigenome.name as epigenome_name from regulatory_build join regulatory_build_epigenome using (regulatory_build_id) join epigenome using (epigenome_id) where regulatory_build.is_current=1',
            },
            -flow_into   => {
               2 => { 'export_regulatory_activities', INPUT_PLUS() },
            },
        },
        {   -logic_name  => 'export_regulatory_activities',
	    -analysis_capacity => 20,
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'export_regulatory_features.pl --epigenome_name "#epigenome_name#" --output_dir #ftp_base_dir#/#species#/RegulatoryFeatureActivity --registry #reg_conf# --species #species#',
            },
        },
        {   -logic_name  => 'export_regulatory_features',
	    -analysis_capacity => 20,
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'export_regulatory_features.pl --only_summary --output_dir #ftp_base_dir#/#species# --registry #reg_conf# --species #species#',
            },
        },
    ]
}

1;
