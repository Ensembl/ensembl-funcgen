package Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::RunSwitchTableEngineHealthchecks_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::Base');

sub pipeline_analyses {
    my $self = shift;
    return [
      {
          -logic_name  => 'start_switch_table_engine_healthchecks',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -flow_into => {
              MAIN => 'hc_switch_table_engine'
          },
      },
      {
          -logic_name  => 'hc_switch_table_engine',
          -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SqlHealthcheck',
          -parameters => {
            db_conn       => 'funcgen:#species#',
            description   => 'Check the tables really have been switched to MyISAM. The switch can fail, if there were foreign key constraints.',
            query         => '
              select 
                table_name, table_rows, update_time, create_time, engine
              from
                information_schema.tables
              where 
                table_schema = database()
                and table_name in (
                  "array",
                  "array_chip",
                  "probe",
                  "probe_feature",
                  "probe_seq",
                  "probe_set",
                  "unmapped_object",
                  "unmapped_reason"
                )
                and engine!="MyISAM"
              order by 
                table_name;
            ',
            expected_size => '0'
          },
      },
    ];
}

1;
