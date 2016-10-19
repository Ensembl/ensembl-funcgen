package Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Funcgen::HiveConfig::Probe2Transcript_conf');

sub _pipeline_analyses_probe_align {
    my $self = shift;
    
    return [
        {   -logic_name  => 'PrePipelineChecks',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::PrePipelineChecks',
            -meadow_type => 'LOCAL',
            -input_ids => [ {} ],
        },
        {
	    -logic_name  => 'MkTmpDir',
            -meadow_type => 'LOCAL',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd => 'mkdir -p #directory#'
            },            
            -input_ids => [ 
	      {
		directory => $self->o('tempdir'),
	      } 
            ],
        },
        {   -logic_name  => 'CreateDB',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::CreateDB',
            -meadow_type => 'LOCAL',
            -input_ids => [ {} ],
            -wait_for => [ 'PrePipelineChecks' ],
        },
        {   -logic_name  => 'JobFactoryImportArrays',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::JobFactory',
            -meadow_type => 'LOCAL',
            # The jobs can only be started after the databases have been created.
            -wait_for => [ 'CreateDB' ],
            -input_ids => [ 
	      {
		probe_directories => $self->o('probe_directories'),
	      },
            ],
	    -flow_into => {
		'1' => 'RollbackArray',
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
               1 => [ 'ImportArrays' ],
            },
        },
        {   -logic_name  => 'ImportArrays',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::ImportArrays',
            -meadow_type => 'LSF',
            -flow_into => {
	       # MEMLIMIT
               -1 => [ 'ImportArrays8Gb' ],
            },
        },
        {   -logic_name  => 'ImportArrays8Gb',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::ImportArrays',
            -meadow_type => 'LSF',
            -rc_name    => '8Gb_job',
            -can_be_empty => 1,
            -flow_into => {
	       # MEMLIMIT
               -1 => [ 'ImportArrays16Gb' ],
            },
            -can_be_empty => 1,
        },
        {   -logic_name  => 'ImportArrays16Gb',
            -meadow_type => 'LSF',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::ImportArrays',
            -rc_name    => '16Gb_job',
            -can_be_empty => 1,
            -flow_into => {
	       # MEMLIMIT
               -1 => [ 'ImportArrays64Gb' ],
            },
            -can_be_empty => 1,
        },
        {   -logic_name  => 'ImportArrays64Gb',
            -meadow_type => 'LSF',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::ImportArrays',
            -rc_name    => '64Gb_job',
            -can_be_empty => 1,
        },
        {   -logic_name  => 'CreateTemporaryIndices',
            -meadow_type => 'LOCAL',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
            # These jobs may fail, if index already exists for some reason.
            -failed_job_tolerance => 100,
            -parameters  => {
              db_conn => $self->_create_db_url_from_dba_hash($self->o('tracking_dba_hash')),
            },
            -input_ids   => [ 
              { sql     => [ 'create index temp_probe_probe_seq_id on probe (probe_seq_id)', ], },
              { sql     => [ 'create index temp_seq_region_probe_analysis_idx on probe_feature (`seq_region_id`,`seq_region_start`, `seq_region_end`, `probe_id`, `analysis_id`)', ], },
            ],
            -wait_for    => [ 'ImportArrays', 'ImportArrays8Gb', 'ImportArrays16Gb', 'ImportArrays64Gb', ],
        },
#         {   -logic_name  => 'DeleteOrphanTrackingData',
#             -meadow_type => 'LOCAL',
#             -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
#             -input_ids   => [ {
#               db_conn => $self->_create_db_url_from_dba_hash($self->o('tracking_dba_hash')),
#               sql     => [ 
#                 # Takes too much time and then blocks table in the meantime
# #                 'delete from probe_seq where probe_seq_id not in (select probe_seq_id from probe)',
# #                 'delete from probe_alias where probe_id not in (select probe_id from probe)',
#                 #
#                 #'delete from probe_set where probe_set_id not in (select distinct probe_set_id from probe)'
#               ],
#             } ],
#             -wait_for    => [ 'CreateTemporaryIndices', ],
#         },
        {   -logic_name  => 'UseInnoDB',
            -meadow_type => 'LOCAL',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
            -input_ids   => [ 
	      {
		db_conn => $self->_create_db_url_from_dba_hash($self->o('tracking_dba_hash')),
		sql     => [ 
		  'ALTER TABLE probe_feature ENGINE=InnoDB;',
		],            
	      } 
            ],
            -wait_for    => [ 'CreateDB', ],
        },        
        {   -logic_name  => 'InsertAnalyses',
            -meadow_type => 'LOCAL',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::InsertAnalyses',
            -wait_for    => [ 'ImportArrays', 'ImportArrays8Gb', 'ImportArrays16Gb', 'ImportArrays64Gb', ],
            -input_ids   => [ {} ],
        },
        {   -logic_name  => 'InsertExternalDb',
            -meadow_type => 'LOCAL',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::InsertExternalDb',
            -wait_for    => [ 'CreateDB', ],
            -input_ids   => [ {} ],
        },
        {   -logic_name  => 'InsertProbeAlias',
            -meadow_type => 'LSF',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::InsertProbeAlias',
            -wait_for => [ 'ImportArrays', 'ImportArrays8Gb', 'ImportArrays16Gb', 'ImportArrays64Gb', 'CreateTemporaryIndices' ],
            -input_ids   => [ {} ],
            -rc_name     => '2Gb_job',
        },
        {   -logic_name => 'DumpToplevel',
            -meadow_type => 'LSF',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd' => 'sequence_dump.pl -dbuser #dbuser# -dbname #dbname# -dbhost #dbhost# -dbport #dbport# -toplevel -onefile -filename #filename# -mask_repeat Dust -mask_repeat RepeatMask',
            },
            -input_ids => [ 
	      {
		filename => $self->o('toplevel_dump_file'),
		dbname   => $self->o('dnadb_name'),
		dbhost   => $self->o('dnadb_host'),
		dbport   => $self->o('dnadb_port'),
		dbuser   => $self->o('dnadb_user'),
	      },
            ],
            -wait_for => [ 'PrePipelineChecks' ],
            -can_be_empty => 1,
            -flow_into => {
	       # MEMLIMIT
               -1 => [ 'DumpToplevel8Gb' ],
            },
        },
        {   -logic_name  => 'DumpToplevel8Gb',
            -meadow_type => 'LSF',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -rc_name     => '8Gb_job',
            -parameters  => {
                'cmd' => 'sequence_dump.pl -dbuser #dbuser# -dbname #dbname# -dbhost #dbhost# -dbport #dbport# -toplevel -onefile -filename #filename# -mask_repeat Dust -mask_repeat RepeatMask',
            },
            -can_be_empty => 1,
        },
        {   -logic_name => 'DumpGenes',
            -meadow_type => 'LSF',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                'cmd' => 'dump_genes.pl -dbuser #dbuser# -dbname #dbname# -dbhost #dbhost# -dbport #dbport# -file #filename# -cdna -stable_id',
            },
            -input_ids => [ 
	      {
		filename => $self->o('transcript_dump_file'),
		dbname   => $self->o('dnadb_name'),
		dbhost   => $self->o('dnadb_host'),
		dbport   => $self->o('dnadb_port'),
		dbuser   => $self->o('dnadb_user'),
	      },
            ],
            -wait_for => [ 'PrePipelineChecks' ],
            -can_be_empty => 1,
            -flow_into => {
	       # MEMLIMIT
               -1 => [ 'DumpGenes8Gb' ],
            },
        },
        {   -logic_name  => 'DumpGenes8Gb',
            -meadow_type => 'LSF',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -rc_name     => '8Gb_job',
            -parameters  => {
                'cmd' => 'dump_genes.pl -dbuser #dbuser# -dbname #dbname# -dbhost #dbhost# -dbport #dbport# -file #filename# -cdna -stable_id',
            },
            -can_be_empty => 1,
        },
        {   -logic_name  => 'DumpUnmappedSeqs',
            -meadow_type => 'LOCAL',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::DumpUnmappedSeqs',
            -input_ids => [ 
	      {
		unmapped_sequences_file => $self->o('unmapped_sequences_file'),
	      },
            ],
            -wait_for => [ 'JobFactoryImportArrays', 'ImportArrays', 'ImportArrays8Gb', 'ImportArrays16Gb', 'ImportArrays64Gb', ],
            -can_be_empty => 1,
        },
        {   -logic_name  => 'JobFactoryProbeAlign',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::FastaFactory',
            -meadow_type => 'LOCAL',
            -input_ids => [ 
	      {
		inputfile => $self->o('unmapped_sequences_file'),
		max_chunk_length => 1000,
		output_dir => $self->o('tempdir'),
		output_prefix => 'probe_chunk_',
		output_suffix => '.fasta',
		hash_directories => 1,
	      }
            ],
            -wait_for => [ 
	      'JobFactoryImportArrays', 
	      'ImportArrays', 'ImportArrays8Gb', 'ImportArrays16Gb', 'ImportArrays64Gb',
	      'DumpToplevel', 'DumpToplevel8Gb', 
	      'DumpGenes', 'DumpGenes8Gb',
	      'DumpUnmappedSeqs', 'CreateTemporaryIndices',
            ],
	    -flow_into => {
		'2' => [ 
		  'ProbeAlignGenomic8Gb',
		  'ProbeAlignTranscript', 
		],
	    },
        },
        {   -logic_name  => 'ProbeAlignGenomic8Gb',
            -meadow_type => 'LSF',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::ProbeAlign',
            -rc_name     => '8Gb_job',
            -can_be_empty => 1,
            -priority => 20,
            -analysis_capacity => 20,
            -batch_size => 100,
            -parameters => {
		QUERYSEQS    => '#chunk_name#',
		OUTDB        => $self->o('tracking_dba_hash'),
		TARGETSEQS   => $self->o('toplevel_dump_file'),
		mapping_type => 'genomic', # genomic or transcript
            },
            -flow_into => {
	       # MEMLIMIT
               -1 => [ 'ProbeAlignGenomic64Gb' ],
            },
        },        
        {   -logic_name  => 'ProbeAlignGenomic64Gb',
	    -can_be_empty => 1,
            -meadow_type => 'LSF',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::ProbeAlign',
            -batch_size => 1,
            -parameters => {
		QUERYSEQS    => '#chunk_name#',
		OUTDB        => $self->o('tracking_dba_hash'),
		TARGETSEQS   => $self->o('toplevel_dump_file'),
		mapping_type => 'genomic', # genomic or transcript
            },
            -rc_name     => '64Gb_job',
            -priority => 30,
            -max_retry_count => 1,
        },        
        {   -logic_name  => 'ProbeAlignTranscript',
	    -can_be_empty => 1,
            -meadow_type => 'LSF',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::ProbeAlign',
            -priority => 10,
            -batch_size => 100,
            -parameters => {
		QUERYSEQS    => '#chunk_name#',
		OUTDB        => $self->o('tracking_dba_hash'),
		TARGETSEQS   => $self->o('transcript_dump_file'),
		mapping_type => 'transcript', # genomic or transcript
            },
            -flow_into => {
	       # MEMLIMIT
               -1 => [ 'ProbeAlignTranscript8Gb' ],
               },
            -wait_for => [ 
              'JobFactoryImportArrays', 
              'ImportArrays', 'ImportArrays8Gb', 'ImportArrays16Gb', 'ImportArrays64Gb',
              'DumpToplevel', 'DumpToplevel8Gb', 
              'DumpGenes', 'DumpGenes8Gb',
              'DumpUnmappedSeqs', 'CreateTemporaryIndices',
            ],
        },
        {   -logic_name  => 'ProbeAlignTranscript8Gb',
	    -can_be_empty => 1,
            -meadow_type => 'LSF',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::ProbeAlign',
            -rc_name     => '8Gb_job',
            -priority => 20,
            -batch_size => 10,
            -parameters => {
		QUERYSEQS    => '#chunk_name#',
		OUTDB        => $self->o('tracking_dba_hash'),
		TARGETSEQS   => $self->o('transcript_dump_file'),
		mapping_type => 'transcript', # genomic or transcript
            },
            -flow_into => {
	       # MEMLIMIT
               -1 => [ 'ProbeAlignTranscript64Gb' ],
            },
        },        
        {   -logic_name  => 'ProbeAlignTranscript64Gb',
            -parameters => {
		QUERYSEQS    => '#chunk_name#',
		OUTDB        => $self->o('tracking_dba_hash'),
		TARGETSEQS   => $self->o('transcript_dump_file'),
		mapping_type => 'transcript', # genomic or transcript
            },
	    -can_be_empty => 1,
            -meadow_type => 'LSF',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::ProbeAlign',
            -rc_name     => '64Gb_job',
            -priority => 30,
            -batch_size => 1,
            -max_retry_count => 1,
            },
        {   -logic_name  => 'DropTemporaryIndices',
            -meadow_type => 'LOCAL',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
            -input_ids   => [ {
              db_conn => $self->_create_db_url_from_dba_hash($self->o('tracking_dba_hash')),
              sql     => [ 
                 'drop index temp_probe_probe_seq_id on probe',
                 'drop index temp_seq_region_probe_analysis_idx on probe_feature',
              ],
            } ],
            -wait_for => [
              'JobFactoryProbeAlign', 
              'ProbeAlignGenomic8Gb', 
              'ProbeAlignGenomic64Gb', 
              'InsertProbeAlias',
              'JobFactoryProbeAlign',
              'ProbeAlignTranscript', 
              'ProbeAlignTranscript8Gb', 
              'ProbeAlignTranscript64Gb', 
            ],
        },
        {   -logic_name  => 'DeleteTempDir',
            -meadow_type => 'LSF',
            -max_retry_count => 1,
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd => 'rm -rf #directory#'
            },
            -input_ids => [ 
	      {
		directory => $self->o('tempdir')
	      },
            ],
            # Must wait for any analysis using the temporary directory
	    -wait_for => [
	      'JobFactoryProbeAlign', 
	      'ProbeAlignGenomic8Gb', 
	      'ProbeAlignGenomic64Gb', 
	      'ProbeAlignTranscript', 
	      'ProbeAlignTranscript8Gb', 
	      'ProbeAlignTranscript64Gb', 
	    ],
        },
    ];    
}

sub pipeline_analyses {
    my $self = shift;
    
    # This does not work!!! Due to the nature of the o mechanism.
    #my $run_probe2transcript_only = $self->o('run_probe2transcript_only');
    
    my $run_probe2transcript_only = 0;
    
    my @all_analyses;
    
    my $pipeline_analyses_probe2transcript = $self->_pipeline_analyses_probe2transcript;
    
    # Connecting analyses as seen in the sub pipeline_analyses in 
    # Bio::EnsEMBL::Compara::PipeConfig::CAFE_conf
    #
    if ($run_probe2transcript_only) {
    
      delete $pipeline_analyses_probe2transcript->[0]->{-wait_for};
      $pipeline_analyses_probe2transcript->[1]->{-wait_for} = [ 
	      'UpdateTranscriptXrefs', 
      ];            
      push @all_analyses, @$pipeline_analyses_probe2transcript;
    
    } else {
    
      my $pipeline_analyses_probe_align = $self->_pipeline_analyses_probe_align;
      
      push @all_analyses, @$pipeline_analyses_probe_align;
      push @all_analyses, @$pipeline_analyses_probe2transcript;
      push @all_analyses, {
	    -logic_name  => 'UseMyIsam',
            -meadow_type => 'LOCAL',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
            -input_ids   => [ {
	      db_conn => $self->_create_db_url_from_dba_hash($self->o('tracking_dba_hash')),
	      sql     => [ 
		'ALTER TABLE probe_feature ENGINE=MyIsam;',
	      ],            
            } ],
            -wait_for    => [ 'Probe2Transcript', ],
        };
    }    
    return \@all_analyses,
}

1;
