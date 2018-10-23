package Bio::EnsEMBL::Funcgen::PipeConfig::ProbeMapping::Populate_meta_coord_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Funcgen::PipeConfig::ProbeMapping::Base');

sub pipeline_analyses {
    my $self = shift;
    
    return [
      {
          -logic_name  => 'start_populate_meta_coords',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -flow_into => {
              MAIN => 'populate_meta_coords'
          },
      },
      {
          -logic_name  => 'populate_meta_coords',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
            cmd => qq( populate_meta_coord.pl    )
              . qq( --species  #species#         )
              . qq( --registry #reg_conf#        )
          },
      },
    ];
}

1;