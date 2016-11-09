package Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::ImportArrays_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Funcgen::HiveConfig::Probe2Transcript_conf');

sub pipeline_analyses {
    my $self = shift;
    
    return [
        {   -logic_name  => 'PrePipelineChecks',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::PrePipelineChecks',
            -input_ids => [ {} ],
            -flow_into => {
                MAIN => 'MkTmpDir',
            },
        },
        {
	    -logic_name  => 'MkTmpDir',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd       => 'mkdir -p #tempdir#',
            },
            -flow_into => {
                MAIN => 'CreateDB',
            },
        },
        {   -logic_name  => 'CreateDB',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::CreateDB',
            -flow_into => {
                MAIN => 'JobFactoryImportArrays',
            },
        },
        {   -logic_name  => 'JobFactoryImportArrays',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::JobFactory',
            -parameters => {
                probe_directories => $self->o('probe_directories'),
            },
	    -flow_into => {
		MAIN => 'RollbackArray',
	    },
        },
        {
	    -logic_name  => 'RollbackArray',
            -meadow_type => 'LOCAL',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -analysis_capacity => 1,
            -parameters => {
                cmd => 'rollback_array.pl'
		  . ' --species '     . $self->o('species')		  
		  . ' -dbhost '       . $self->o('tracking_host')
		  . ' -dbname '       . $self->o('tracking_dbname')
		  . ' -dbuser '       . $self->o('tracking_user')
		  . ' -dbport '       . $self->o('tracking_port')
		  
		  . ( $self->o('tracking_pass') ? ' -dbpass ' . $self->o('tracking_pass') : '' )
		  
		  . ' -arrays #all_array_names# -force',
            },
            -flow_into => {
               '1->A' => 'ImportArrays',
               'A->1' => 'hc_import_arrays',
            },
        },
        {   -logic_name  => 'hc_import_arrays',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SqlHealthcheck',
            -parameters  => {
              query         => 'select distinct array.name from array join array_chip using (array_id) left join probe using (array_chip_id) where probe.probe_id is null and array.class="#array_format#"',
              expected_size => '=0',
              db_conn       => $self->_create_db_url_from_dba_hash($self->o('tracking_dba_hash')),
              description   => 'Check that the arrays imported from #array_format# have probes in the database.';
            },
        },
        {   -logic_name  => 'ImportArrays',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::ImportArrays',
            -flow_into => {
               # MEMLIMIT
               -1 => [ 'ImportArrays8Gb' ],
            },
        },
        {   -logic_name  => 'ImportArrays8Gb',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::ImportArrays',
            -rc_name    => '8Gb_job',
            -flow_into => {
	       # MEMLIMIT
               -1 => [ 'ImportArrays16Gb' ],
            },
        },
        {   -logic_name  => 'ImportArrays16Gb',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::ImportArrays',
            -rc_name    => '16Gb_job',
            -flow_into => {
	       # MEMLIMIT
               -1 => [ 'ImportArrays64Gb' ],
            },
        },
        {   -logic_name  => 'ImportArrays64Gb',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::ImportArrays',
            -rc_name    => '64Gb_job',
        },
    ];
}

1;
