package Bio::EnsEMBL::Funcgen::Ftp::PipeConfig::ExportRegulatoryFeatures;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Hive::PipeConfig::EnsemblGeneric_conf';

sub default_options {
    my ($self) = @_;
    return {
        %{ $self->SUPER::default_options() },
	
	ontology_database => {
	  -species => 'multi',
	  -group   => 'ontology',
	  -host    => 'ens-staging2',
	  -user    => 'ensro',
	  -pass    => undef,
	  -driver  => 'mysql',
	  -port    => 3306,
	  -dbname  => 'ensembl_ontology_85',
	  -db_version => $self->o('ensembl_release')
	},
	regulation_database => {
	  -species => 'homo_sapiens',
	  -group   => 'funcgen',
	  -host    => 'ens-genomics2',
	  -user    => 'ensro',
	  -pass    => undef,
	  -driver  => 'mysql',
	  -port    => 3306,
	  -dbname  => 'mn1_e_homo_sapiens_funcgen_85_38',
	  -db_version => $self->o('ensembl_release')
	},
    };
}

sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
      %{$self->SUPER::pipeline_wide_parameters},

      regulation_database => $self->o('regulation_database'),
      ontology_database   => $self->o('ontology_database'),

      regulation_database_url => $self->dbconn_2_url('regulation_database'),
      ontology_database_url   => $self->dbconn_2_url('ontology_database'),

      ftp_base_dir        => '/lustre/scratch109/ensembl/funcgen/mn1/ersa/mn1_dev3_homo_sapiens_funcgen_85_38/output/mn1_dev3_homo_sapiens_funcgen_85_38/ftp_site',
    };
}

sub pipeline_analyses {
    my $self = shift;

    return [
        {   -logic_name  => 'job_factory',
            -meadow_type => 'LOCAL',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -parameters  => { db_conn => '#regulation_database#' },
            -input_ids   => [ {inputquery => 'select name as epigenome_name from epigenome' } ],
            -flow_into   => {
               2 => 'export_regulatory_features',
            },
        },
        {   -logic_name  => 'export_annotated_features',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -input_ids   => [ {} ],
            -parameters  => {
                cmd => 'export_annotated_features.pl --ftp_base_dir #ftp_base_dir# --regulation_database_url #regulation_database_url# --ontology_database_url #ontology_database_url#',

            },
        },
        {   -logic_name  => 'export_motif_features',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -input_ids   => [ {} ],
            -parameters  => {
                cmd => 'export_motif_features.pl --ftp_base_dir #ftp_base_dir# --regulation_database_url #regulation_database_url# --ontology_database_url #ontology_database_url#',

            },
        },
        {   -logic_name  => 'export_regulatory_features',
	    -analysis_capacity => 20,
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'export_regulatory_features.pl --epigenome_name "#epigenome_name#" --ftp_base_dir #ftp_base_dir# --regulation_database_url #regulation_database_url# --ontology_database_url #ontology_database_url#',
            },
        },
    ]
}

1;
