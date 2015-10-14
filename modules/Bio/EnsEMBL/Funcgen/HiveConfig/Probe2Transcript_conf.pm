package Bio::EnsEMBL::Funcgen::HiveConfig::Probe2Transcript_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::PipeConfig::EnsemblGeneric_conf');

use Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::Utils qw (create_db_url_from_dba_hash);

sub _hash_placeholder {
  return { 'No final value available yet.' => 'No final value available yet.' };
}

sub _create_db_url_from_dba_hash {

  my $self = shift;  
  my $hash = shift;  
  if (ref $hash eq 'HASH') {
    return create_db_url_from_dba_hash($hash);
  }
  # Allows dereferencing
  return _hash_placeholder;
}

sub _dereference_hash {
  my $hashref = shift;  
  if (ref $hashref eq 'HASH') {
    return %$hashref;
  }
  return %{&_hash_placeholder};
}

sub beekeeper_extra_cmdline_options {
    my ($self) = @_;
    return '-keep_alive -can_respecialize 1';
}

sub default_options {
    my ($self) = @_;
    return {
        %{ $self->SUPER::default_options() },

	dnadb_name           => 'homo_sapiens_core_81_38',
	dnadb_host           => 'ens-livemirror',
	dnadb_port           => '3306',
	dnadb_user           => 'ensro',
	dnadb_pass           => '',

	tracking_user        => 'ensadmin',
	tracking_pass        => 'ensembl',
	tracking_dbname      => 'mn1_tracking_homo_sapiens_funcgen_81_3',
	tracking_port        => '3306',
	tracking_host        => 'ens-genomics2',
	
	species              => 'homo_sapiens',
	
	probe_directories    => '/lustre/scratch110/ensembl/funcgen/array_mapping/HOMO_SAPIENS/',

	toplevel_dump_file      => $self->o('tempdir') . '/target.fasta',
	transcript_dump_file    => $self->o('tempdir') . '/target_genes.fasta',
	unmapped_sequences_file => $self->o('tempdir') . '/unmapped_probe_sequences.fasta',
	
	tempdir => '/lustre/scratch110/ensembl/funcgen/array_mapping/'.$ENV{USER}.'/temp/' .$self->o('species'),
	
	tracking_dba_hash => {
	    -user         => $self->o('tracking_user'),
	    -pass         => $self->o('tracking_pass'),
	    -port         => $self->o('tracking_port'),
	    -dbname       => $self->o('tracking_dbname'),
	    -host         => $self->o('tracking_host'),

	    -dnadb_name   => $self->o('dnadb_name'),
	    -dnadb_host   => $self->o('dnadb_host'),
	    -dnadb_port   => $self->o('dnadb_port'),
	    -dnadb_user   => $self->o('dnadb_user'),
	    -dnadb_pass   => $self->o('dnadb_pass'),

	    -species      => $self->o('species'),
	},
    };
}

sub pipeline_wide_parameters {
    my ($self) = @_;
	
    return {
        %{$self->SUPER::pipeline_wide_parameters},     
 	pipeline_repository_dir => $self->o('pipeline_repository_dir'),
	tracking_dba_hash => $self->o('tracking_dba_hash'),
	tempdir => $self->o('tempdir'),
    };
}

sub _pipeline_analyses_probe2transcript {
    my $self = shift;

    return [
        {   -logic_name  => 'UpdateTranscriptXrefs',
            -meadow_type => 'LSF',
            -max_retry_count => 1,
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -can_be_empty => 1,
            -rc_name    => '4Gb_job',
            -parameters => {
                'cmd' => 'update_transcript_xrefs.pl --species #species# --transcript_dbname #dnadb_name# --transcript_host #dnadb_host# --transcript_port #dnadb_port# --transcript_user #dnadb_user# #dnadb_pass# --xref_host #xref_host# --xref_dbname #xref_dbname# --xref_user #xref_user# #xref_pass# --xref_port #xref_port#',
            },
            -input_ids => [ 
	      {           
		filename => $self->o('transcript_dump_file'),
		species => $self->o('species'),

		xref_dbname   => $self->o('tracking_dbname'),
		xref_host     => $self->o('tracking_host'),
		xref_port     => $self->o('tracking_port'),
		xref_user     => $self->o('tracking_user'),
		xref_pass     => ($self->o('tracking_pass') ? '--xref_pass ' . $self->o('tracking_pass') : '' ),

		dnadb_name     => $self->o('dnadb_name'),
		dnadb_host     => $self->o('dnadb_host'),
		dnadb_port     => $self->o('dnadb_port'),
		dnadb_user     => $self->o('dnadb_user'),		
		dnadb_pass     => ($self->o('dnadb_pass') ? '--transcript_pass ' . $self->o('tracking_pass') : '' ),
	      },
            ],
            -wait_for => [ 
	      'InsertExternalDb', 
	      
	      # UpdateTranscriptXrefs needs a lot of exclusive access to the 
	      # database tables. This waiting rule prevents it from running
	      # at the beginning, when other analyses are also heavily 
	      # accessing the database.
	      #
	      'JobFactoryProbeAlign',
            ],            
        },
        {   -logic_name  => 'P2TJobFactory',
            -meadow_type => 'LOCAL',
            -max_retry_count => 1,
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -can_be_empty => 1,
            -parameters => {
                db_conn => $self->_create_db_url_from_dba_hash($self->o('tracking_dba_hash')),
            },
            -input_ids => [ 
	      {
 		inputquery => 'select group_concat(name separator " ") as arrays, vendor, class from array where format!="METHYLATION" group by vendor, class',
	      },
            ],
            -flow_into => {
               2 => [ 'Probe2Transcript' ],
            },
            -wait_for => [ 
	      'UpdateTranscriptXrefs', 
	      'ProbeAlignTranscript', 
	      'ProbeAlignTranscript8Gb', 
	      'ProbeAlignTranscript64Gb', 
            ],
        },
        {   -logic_name  => 'Probe2Transcript',
            -meadow_type => 'LSF',
            -max_retry_count => 1,
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -rc_name    => '16Gb_job',
            -parameters => {
            
                cmd => 'probe2transcript.pl'

		  . ' --calculate_utrs'
		  . ' --utr_multiplier 1'
		  
		  . ' --filename ' . 'probe2transcript.' . $self->o('species') . '.#vendor#.log'

		  . ' --species '           . $self->o('species')
		  
		  . ' --transcript_dbname ' . $self->o('dnadb_name')
		  . ' --transcript_host '   . $self->o('dnadb_host')
		  . ' --transcript_port '   . $self->o('dnadb_port')
		  . ' --transcript_user '   . $self->o('dnadb_user')
		  . ( $self->o('dnadb_pass') ? ' --transcript_pass ' . $self->o('dnadb_pass') : '')
		  
		  . ' --xref_host '         . $self->o('tracking_host')
		  . ' --xref_dbname '       . $self->o('tracking_dbname')
		  . ' --xref_user '         . $self->o('tracking_user')
		  . ' --xref_port '         . $self->o('tracking_port')
		  . ( $self->o('tracking_pass') ? ' --xref_pass ' . $self->o('tracking_pass') : '' )
		  
		  . ' --arrays #arrays# -vendor #vendor# -format #class#',
            },
	    -wait_for => [
	    
	      # Not really a dependency, but UpdateTranscriptXrefs puts a lot of 
	      # load on the server, so it is best to wait until it is over.
	      #
	      'UpdateTranscriptXrefs', 
	    ],
        },
    ]
}

sub pipeline_analyses {
    my $self = shift;
    
    my $pipeline_analyses_probe2transcript = $self->_pipeline_analyses_probe2transcript;

    my @all_analyses;
    delete $pipeline_analyses_probe2transcript->[0]->{-wait_for};
    $pipeline_analyses_probe2transcript->[1]->{-wait_for} = [ 
	    'UpdateTranscriptXrefs', 
    ];            
    push @all_analyses, @$pipeline_analyses_probe2transcript;
        
    return \@all_analyses,
}

# Copied over from Bio::EnsEMBL::Compara::PipeConfig::ProteinTrees_conf
#
# Added -q long to avoid RUNLIMIT problems for 64Gb jobs
#
sub resource_classes {
    my ($self) = @_;
    return {
        %{$self->SUPER::resource_classes},  # inherit 'default' from the parent class

         '250Mb_job'    => {'LSF' => '-C0 -M250   -R"select[mem>250]   rusage[mem=250]"' },
         '500Mb_job'    => {'LSF' => '-C0 -M500   -R"select[mem>500]   rusage[mem=500]"' },
         '1Gb_job'      => {'LSF' => '-C0 -M1000  -R"select[mem>1000]  rusage[mem=1000]"' },
         '2Gb_job'      => {'LSF' => '-C0 -M2000  -R"select[mem>2000]  rusage[mem=2000]"' },
         '4Gb_job'      => {'LSF' => '-C0 -M4000  -R"select[mem>4000]  rusage[mem=4000]"' },
         '8Gb_job'      => {'LSF' => '-C0 -M8000  -R"select[mem>8000]  rusage[mem=8000]"' },
         '16Gb_job'     => {'LSF' => '-C0 -M16000 -R"select[mem>16000] rusage[mem=16000]"' },
         '24Gb_job'     => {'LSF' => '-C0 -M24000 -R"select[mem>24000] rusage[mem=24000]"' },
         '32Gb_job'     => {'LSF' => '-C0 -M32000 -R"select[mem>32000] rusage[mem=32000]"' },
         '48Gb_job'     => {'LSF' => '-C0 -M48000 -R"select[mem>48000] rusage[mem=48000]"' },
         '64Gb_job'     => {'LSF' => '-C0 -M64000 -R"select[mem>64000] rusage[mem=64000]" -q long' },

         '16Gb_16c_job' => {'LSF' => '-n 16 -C0 -M16000 -R"select[mem>16000] rusage[mem=16000]"' },
         '64Gb_16c_job' => {'LSF' => '-n 16 -C0 -M64000 -R"select[mem>64000] rusage[mem=64000]"' },

         '8Gb_64c_mpi'  => {'LSF' => '-q mpi -n 64 -a openmpi -M8000 -R"select[mem>8000] rusage[mem=8000] same[model] span[ptile=16]"' },
         '32Gb_64c_mpi' => {'LSF' => '-q mpi -n 64 -a openmpi -M32000 -R"select[mem>32000] rusage[mem=32000] same[model] span[ptile=16]"' },

    };
}

1;



