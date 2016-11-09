package Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::ImportArraysNew_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::Base');

=head1

  export HIVE_URL='mysql://ensadmin:ensembl@ens-genomics2:3306/mmn1_tracking_homo_sapiens_funcgen_87_38_hive'
  ftp_pipeline_parameters="-pipeline_url $HIVE_URL"
  
  init_pipeline.pl Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::ImportArraysNew_conf $ftp_pipeline_parameters -hive_force_init 1

=cut

sub pipeline_wide_parameters {
    my $self = shift;
    return {
      %{$self->SUPER::pipeline_wide_parameters},

      tempdir          => $self->o('tempdir'),
      probe_directory  => $self->o('probe_directory'),
      reg_conf         => $self->o('reg_conf'),
    };
}

sub default_options {
    my ($self) = @_;
    return {
      %{ $self->SUPER::default_options() },

      tempdir          => '/lustre/scratch109/ensembl/funcgen/array_mapping_temp',
      probe_directory  => '/lustre/scratch109/ensembl/funcgen/array_mapping/',
      reg_conf         => '/nfs/users/nfs_m/mn1/work_dir_probemapping/lib/ensembl-funcgen/registry.pm',
    };
}

sub beekeeper_extra_cmdline_options {
    my ($self) = @_;
    return '-reg_conf ' . $self->o('reg_conf') . ' -keep_alive -can_respecialize 1';
}

sub pipeline_analyses {
    my $self = shift;
    
    return [
        {
            -logic_name  => 'start',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -input_ids => [
              { species => 'homo_sapiens', },
            ],
            -flow_into => {
                MAIN => 'pre_pipeline_checks',
            },
        },
        {   -logic_name  => 'pre_pipeline_checks',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::PrePipelineChecks',
            -flow_into => {
                MAIN => 'make_temp_dir',
            },
        },
        {
            -logic_name  => 'make_temp_dir',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd       => 'mkdir -p #tempdir#/#species#',
            },
            -flow_into => {
                MAIN => 'rollback_array',
            },
        },
        {
            -logic_name  => 'rollback_array',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
            -parameters => {
                sql     => [
                  "truncate array;",
                  "truncate array_chip;",
                  "truncate probe;",
                  "truncate probe_feature;",
                  "truncate probe_alias;",
                  "truncate probe_seq;",
                  "truncate probe_set;",
                ],
                db_conn => 'funcgen:#species#',
            },
            -flow_into => {
               MAIN => 'job_factory_import_arrays',
            },
        },
        {
          -logic_name  => 'job_factory_import_arrays',
          -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::JobFactory',
          -parameters => {
              probe_directories => '#probe_directory#/#species#',
          },
          -flow_into => {
            MAIN => 'parse_probe_fasta_file',
          },
        },
        {
            -logic_name  => 'parse_probe_fasta_file',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd       => '
                  import_parse_probe_fasta_file.pl \
                    --array_name      #array_class# \
                    --probe_file      #probe_file# \
                    --parsed_output   #tempdir#/#species#/#array_class#_parsed_probes.pl
                ',
            },
            -flow_into => {
                MAIN => 'create_array_objects',
            },
        },
        {
            -logic_name  => 'create_array_objects',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd       => '
                  import_create_array_objects.pl \
                    --array_name        #array_class# \
                    --parsed_probe_data #tempdir#/#species#/#array_class#_parsed_probes.pl \
                    --output_file       #tempdir#/#species#/#array_class#_array_objects.pl
                  ',
            },
            -flow_into => {
                MAIN => 'store_array_objects',
            },
        },
        {
            -logic_name  => 'store_array_objects',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -analysis_capacity => 1,
            -parameters => {
                cmd       => '
                  import_store_array_objects.pl \
                    --registry           #reg_conf# \
                    --species            #species# \
                    --array_objects_file #tempdir#/#species#/#array_class#_array_objects.pl
                ',
            },
          -flow_into => {
              MEMLIMIT => 'store_array_objects_himem',
              MAIN     => 'insert_probe_alias',
          },
        },
        {
            -logic_name  => 'store_array_objects_himem',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -analysis_capacity => 1,
            -parameters => {
                cmd       => '
                  import_store_array_objects.pl \
                    --registry           #reg_conf# \
                    --species            #species# \
                    --array_objects_file #tempdir#/#species#/#array_class#_array_objects.pl
                ',
            },
            -rc_name     => '16Gb_job',
            -flow_into => {
                MAIN     => 'insert_probe_alias',
            },
        },
        {
            -logic_name  => 'insert_probe_alias',
            -meadow_type => 'LSF',
            -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::InsertProbeAlias',
            -rc_name     => '2Gb_job',
        },
    ];
}

sub resource_classes {
    my ($self) = @_;
    return {
        %{$self->SUPER::resource_classes},  # inherit 'default' from the parent class

        'default' => {
          'LSF'   => ['', '--reg_conf '.$self->o('reg_conf')], 
          'LOCAL' => ['', '--reg_conf '.$self->o('reg_conf')] 
        },
        '250Mb_job'    => {'LSF' => '-M250   -R"select[mem>250]   rusage[mem=250]"' },
        '500Mb_job'    => {'LSF' => '-M500   -R"select[mem>500]   rusage[mem=500]"' },
        '1Gb_job'      => {'LSF' => '-M1000  -R"select[mem>1000]  rusage[mem=1000]"' },
        '2Gb_job'      => {'LSF' => '-M2000  -R"select[mem>2000]  rusage[mem=2000]"' },
        '4Gb_job'      => {'LSF' => '-M4000  -R"select[mem>4000]  rusage[mem=4000]"' },
        '8Gb_job'      => {'LSF' => '-M8000  -R"select[mem>8000]  rusage[mem=8000]"' },
        '16Gb_job'     => {
          'LSF' => [ 
            '-M16000 -R"select[mem>16000] rusage[mem=16000]"', 
            '--reg_conf ' . $self->o('reg_conf') 
          ]
        },
        '24Gb_job'     => {'LSF' => '-M24000 -R"select[mem>24000] rusage[mem=24000]"' },
        '32Gb_job'     => {'LSF' => '-M32000 -R"select[mem>32000] rusage[mem=32000]"' },
        '48Gb_job'     => {'LSF' => '-M48000 -R"select[mem>48000] rusage[mem=48000]"' },
        '64Gb_job'     => {'LSF' => '-M64000 -R"select[mem>64000] rusage[mem=64000]"' },
    };
}


1;
