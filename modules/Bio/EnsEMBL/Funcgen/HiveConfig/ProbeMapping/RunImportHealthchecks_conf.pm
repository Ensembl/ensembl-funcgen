package Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::RunImportHealthchecks_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::Base');

sub pipeline_analyses {
    my $self = shift;
    
    return [
      {
          -logic_name  => 'start_import_healthchecks',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -flow_into => {
              MAIN => 'hc_probes_link_to_sequence',
          },
      },

      {
          -logic_name  => 'hc_probes_link_to_sequence',
          -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SqlHealthcheck',
          -parameters => {
            db_conn       => 'funcgen:#species#',
            description   => 'Assert all probes have a probe sequence',
            query         => 'select probe.probe_id, array.class from probe join array_chip using (array_chip_id) join array using (array_id) where probe_seq_id is null',
            expected_size => '0'
          },
        -flow_into => {
            MAIN => 'hc_probe_sequence_not_empty',
        },
      },
      {
          -logic_name  => 'hc_probe_sequence_not_empty',
          -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SqlHealthcheck',
          -parameters => {
            db_conn       => 'funcgen:#species#',
            description   => 'Assert probe sequences are not empty strings',
            query         => 'select probe_seq.sequence, array.class from probe_seq join probe using (probe_seq_id) join array_chip using (array_chip_id) join array using (array_id) where sequence = ""',
            expected_size => '0'
          },
       -flow_into => {
           MAIN => 'hc_probe_set_sizes_ok',
       },
      },
      {
          -logic_name  => 'hc_probe_set_sizes_ok',
          -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SqlHealthcheck',
          -parameters => {
            db_conn       => 'funcgen:#species#',
            description   => 'Assert probe sets are set to the correct size',
            query         => '
                select * from probe_set join (
                        select probe_set_id, array_chip_id, count(*) as counted_size from probe group by probe_set_id, array_chip_id
                ) as count_them using (probe_set_id, array_chip_id) where count_them.counted_size != probe_set.size
            ',
            expected_size => '0'
          },
      },
    ];
}

1;
