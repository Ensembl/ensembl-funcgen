package Bio::EnsEMBL::Funcgen::HiveConfig::Ftp::ExportAnnotatedFeatures;

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
	  -species => 'mus_musculus',
	  -group   => 'funcgen',
	  -host    => 'ens-staging2',
	  -user    => 'ensro',
	  -pass    => undef,
	  -driver  => 'mysql',
	  -port    => 3306,
	  -dbname  => 'mus_musculus_funcgen_85_38',
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

      ftp_base_dir  => '/lustre/scratch109/ensembl/funcgen/mn1/ersa/mn1_dev3_homo_sapiens_funcgen_85_38/output/mn1_dev3_homo_sapiens_funcgen_85_38/ftp_site/mus_musculus',
      temp_dir      => '#ftp_base_dir#/tempdir/annotated_features',
    };
}

sub pipeline_analyses {
    my $self = shift;

    return [
        {   -logic_name  => 'job_factory_annotated_features',
            -module      => 'Bio::EnsEMBL::Funcgen::Hive::Ftp::JobFactoryAnnotatedFeatures',
            -parameters  => { 
	      db_conn  => '#regulation_database#',
            },
            -input_ids   => [ {} ],
            -flow_into   => {
               '2->A' => 'export_annotated_features',
               'A->1' => 'merge_annotated_features',
            },
        },
        {   -logic_name  => 'export_annotated_features',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -analysis_capacity => 10,
            -batch_size => 50,
            -parameters  => {
                cmd => 'export_annotated_features.pl --output_file #directory#/#file# --regulation_database_url #regulation_database_url# --ontology_database_url #ontology_database_url# --min_id #min_id# --max_id #max_id#',
            },
        },
        {   -logic_name  => 'merge_annotated_features',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'cat #gff_files_from_batches# | xargs --max-args=50 cat >> #merged_gff#',
            },
            -flow_into   => {
               MAIN => 'gzip_annotated_features',
            },
        },
        {   -logic_name  => 'gzip_annotated_features',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'gzip #merged_gff#',
            },
            -flow_into   => {
               MAIN => 'mv_annotated_features_to_ftp',
            },
        },
        {   -logic_name  => 'mv_annotated_features_to_ftp',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'mv #merged_gff#.gz #ftp_base_dir#',
            },
            -flow_into   => {
               MAIN => 'rm_annotated_features_temp_dir',
            },
        },
        {   -logic_name  => 'rm_annotated_features_temp_dir',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'rm -rf #temp_dir#',
            },
        },
    ]
}

1;
