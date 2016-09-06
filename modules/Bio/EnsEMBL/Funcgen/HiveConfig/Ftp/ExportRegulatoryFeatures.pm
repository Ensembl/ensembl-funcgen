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
        {   -logic_name  => 'split_species_string',
            -module      => 'Bio::EnsEMBL::Funcgen::Hive::Ftp::SplitString',
            -parameters  => { 
	      string         => '#species_list#',
	      separator      => ',',
	      list_item_name => 'species'
            },
            -input_ids   => [{
	      species_list => $self->o('species_list')
            }
            ],
            -flow_into   => {
               2 => [ 
		'dbconn_for_species', 
		'export_motif_features' 
               ],
            },
        },
        {   -logic_name  => 'dbconn_for_species',
            -module      => 'Bio::EnsEMBL::Funcgen::Hive::Ftp::DbconnForSpecies',
            -parameters  => { 
	      group => 'funcgen'
            },
            -flow_into   => {
               2 => 'job_factory_regulatory_features',
            },
        },
        {   -logic_name  => 'job_factory_regulatory_features',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -parameters  => { 
	      db_conn    => '#url#',
	      inputquery => 'select name as epigenome_name from epigenome',
            },
            -flow_into   => {
               2 => { 'export_regulatory_features', INPUT_PLUS() },
            },
        },
        {   -logic_name  => 'export_motif_features',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'export_motif_features.pl --ftp_base_dir #ftp_base_dir#/#species# --registry #reg_conf# --species #species#',

            },
        },
        {   -logic_name  => 'export_regulatory_features',
	    -analysis_capacity => 20,
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'export_regulatory_features.pl --epigenome_name "#epigenome_name#" --ftp_base_dir #ftp_base_dir#/#species# --registry #reg_conf# --species #species#',
            },
        },
    ]
}

1;
