package Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::SwitchToMyIsam_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::Base');

sub pipeline_analyses {
    my $self = shift;
    
    return [
      {
          -logic_name  => 'start_switch_table_engines',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -flow_into => {
              MAIN => 'switch_to_myisam'
          },
      },
      {
          -logic_name  => 'switch_to_myisam',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
          -parameters => {
              sql     => [
                "ALTER TABLE array           ENGINE=MyISAM;",
                "ALTER TABLE array_chip      ENGINE=MyISAM;",
                "ALTER TABLE probe           ENGINE=MyISAM;",
                "ALTER TABLE probe_feature   ENGINE=MyISAM;",
                "ALTER TABLE probe_seq       ENGINE=MyISAM;",
                "ALTER TABLE probe_set       ENGINE=MyISAM;",
                "ALTER TABLE unmapped_object ENGINE=MyISAM;",
                "ALTER TABLE unmapped_reason ENGINE=MyISAM;",
              ],
              db_conn => 'funcgen:#species#',
          },
      },
    ];
}

1;
