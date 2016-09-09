package Bio::EnsEMBL::Funcgen::HiveConfig::Ftp::Backbone;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Funcgen::HiveConfig::Ftp::Base';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub pipeline_analyses {
    my $self = shift;

    return [
        {   -logic_name  => 'backbone_fire_exports',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'start_export',
               'A->1' => 'backbone_fire_createmetafiles'
            },
        },
        {   -logic_name  => 'backbone_fire_createmetafiles',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into   => {
               '1->A' => 'create_readme_files',
               'A->1' => 'backbone_pipeline_finished'
            },
        },
        {   -logic_name => 'backbone_pipeline_finished',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        }
    ]
}

1;
