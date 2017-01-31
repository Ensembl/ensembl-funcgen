package Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::ImportArraysNew_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::Base');

=head1

  export HIVE_URL='mysql://ensadmin:ensembl@ens-genomics2:3306/mmn1_tracking_homo_sapiens_funcgen_87_38_hive'
  ftp_pipeline_parameters="-pipeline_url $HIVE_URL"
  
  init_pipeline.pl Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::ImportArraysNew_conf $ftp_pipeline_parameters -hive_force_init 1

=cut

sub pipeline_analyses {
    my $self = shift;
    
    return [
        {
            -logic_name  => 'start_import',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
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
                  "truncate probe_seq;",
                  "truncate probe_set;",
                  "truncate probe_transcript;",
                  "truncate probeset_transcript;",
                  "truncate probe_feature_transcript;",
                  "truncate unmapped_object;",
                  "truncate unmapped_reason;",
                  "delete analysis_description from analysis_description, analysis where analysis.analysis_id=analysis_description.analysis_id and logic_name like '%Probe%Align';",
                  "delete from analysis where logic_name like '%Probe%Align';",
                  "delete analysis_description from analysis_description, analysis where analysis.analysis_id=analysis_description.analysis_id and logic_name = 'probe2transcript';",
                  "delete from analysis where logic_name = 'probe2transcript';",
                  
                ],
                db_conn => 'funcgen:#species#',
            },
            -flow_into => {
               MAIN => 'create_probe_mapping_analyses',
            },
        },
        {
            -logic_name  => 'create_probe_mapping_analyses',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd       => 
                    'create_probe_mapping_analyses.pl'
                  . ' --registry #reg_conf#'
                  . ' --species  #species#'
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
        },
    ];
}

1;
