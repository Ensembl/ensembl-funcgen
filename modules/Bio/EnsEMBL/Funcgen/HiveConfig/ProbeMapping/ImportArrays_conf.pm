package Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::ImportArrays_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::Base');

sub pipeline_analyses {
    my $self = shift;
    
    return [
        {
            -logic_name  => 'start_import',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into => {
                MAIN => 'truncate_array_tables',
            },
        },
        {
            -logic_name  => 'truncate_array_tables',
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
                  "truncate probe_set_transcript;",
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
               MAIN => 'switch_array_tables_to_innodb',
            },
        },
        {
            -logic_name  => 'switch_array_tables_to_innodb',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
            -parameters => {
                sql     => [
                  "ALTER TABLE array           ENGINE=InnoDB;",
                  "ALTER TABLE array_chip      ENGINE=InnoDB;",
                  "ALTER TABLE probe           ENGINE=InnoDB;",
                  "ALTER TABLE probe_feature   ENGINE=InnoDB;",
                  "ALTER TABLE probe_seq       ENGINE=InnoDB;",
                  "ALTER TABLE probe_set       ENGINE=InnoDB;",
                  "ALTER TABLE unmapped_object ENGINE=InnoDB;",
                  "ALTER TABLE unmapped_reason ENGINE=InnoDB;",
                  # Not converting:
                  #
                  # probe_transcript
                  # probe_set_transcript
                  # probe_feature_transcript
                  #
                  # because they are populated by load statements. These are 
                  # unlikely to improve by using innodb.
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
#             MAIN => 'parse_probe_fasta_file',
#             'MAIN->A' => 'parse_probe_fasta_file',
#             'A->MAIN' => 'import_arrays_done',

            '2->A' => 'parse_probe_fasta_file',
            'A->1' => 'import_arrays_done',
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
            -analysis_capacity => 70,
            -parameters => {
                cmd       => '
                  import_store_array_objects.pl \
                    --registry           #reg_conf# \
                    --species            #species# \
                    --array_objects_file #tempdir#/#species#/#array_class#_array_objects.pl
                ',
            },
        },
        {
            -logic_name  => 'import_arrays_done',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into => {
                MAIN => 'run_sql_to_fix_probe_set_issues',
            },
        },

  {
      -logic_name  => 'run_sql_to_fix_probe_set_issues',
      -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
      -analysis_capacity => 1,
      -parameters => {
        db_conn       => 'funcgen:#species#',
        sql           => [
          'drop table if exists probe_set_fixed;',
          '
          create table probe_set_fixed (
            probe_set_id     int(10) unsigned NOT NULL AUTO_INCREMENT,
            probe_set_id_old int(10),
            name             varchar(100) NOT NULL,
            array_chip_id    int(10),
            size             smallint(6) unsigned NOT NULL,
            family           varchar(20) DEFAULT NULL,
            PRIMARY KEY (probe_set_id),
            KEY name (name)
          );
          ',
          '   
          insert into probe_set_fixed (probe_set_id_old, name, array_chip_id, size) 
          select 
            probe_set.probe_set_id as probe_set_id_old, 
            probe_set.name, 
            probe.array_chip_id, 
            count(distinct probe.probe_id) as size
          from 
            probe join probe_set using (probe_set_id) 
            group by 
            probe_set.probe_set_id, 
            probe_set.name, 
            probe.array_chip_id
          ;
          ',
          # Without the index, the next update can be very slow on human
          'create index temp_probe_index on probe (probe_set_id,array_chip_id);',
          '
          update 
            probe, probe_set_fixed 
          set 
            probe.probe_set_id = probe_set_fixed.probe_set_id
          where
            probe.probe_set_id=probe_set_fixed.probe_set_id_old
            and probe.array_chip_id=probe_set_fixed.array_chip_id
          ;
          ',
          'truncate probe_set;',
          '
          insert into probe_set (probe_set_id, name, array_chip_id, size) 
          select 
          probe_set_id, name, array_chip_id, size
          from 
          probe_set_fixed
          ;
          ',
          'drop index temp_probe_index on probe;',
          'drop table probe_set_fixed;'
        ],
      },
  },




    ];
}

1;
