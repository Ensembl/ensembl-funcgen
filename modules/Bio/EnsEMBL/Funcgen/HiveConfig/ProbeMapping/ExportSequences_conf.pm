package Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::ExportSequences_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::Base');

sub pipeline_analyses {
    my $self = shift;
    
    return [
      {
          -logic_name  => 'start_export',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -flow_into => {
              MAIN => [
                'connection_details_as_parameters',
                'export_unmapped_sequences'
              ]
          },
      },
      {
          -logic_name  => 'connection_details_as_parameters',
          -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::ConnectionDetailsAsParameters',
          -parameters => {
              type => 'core',
          },
          -flow_into => {
              MAIN => [
                'export_toplevel',
                'export_transcript_sequences',
              ],
          },
      },
      {   -logic_name => 'export_toplevel',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
              'cmd' => 'sequence_dump.pl '
              . ' -dbuser #username# '
              . ' -dbname #dbname# '
              . ' -dbhost #host# '
              . ' -dbport #port# '
              . ' -toplevel '
              . ' -onefile '
              . ' -filename #toplevel_sequences_file# '
              . ' -mask_repeat Dust '
              . ' -mask_repeat RepeatMask'
              . ' #expr( #password# ? " -dbpass " . #password# : "" )expr# ',
          },
          -rc_name   => '4Gb_job',
          -flow_into => {
              MEMLIMIT => 'export_toplevel_himem',
          },
      },
      {   -logic_name => 'export_toplevel_himem',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -rc_name   => '16Gb_job',
          -parameters => {
              'cmd' => 'sequence_dump.pl '
              . ' -dbuser #username# '
              . ' -dbname #dbname# '
              . ' -dbhost #host# '
              . ' -dbport #port# '
              . ' -toplevel '
              . ' -onefile '
              . ' -filename #toplevel_sequences_file# '
              . ' -mask_repeat Dust '
              . ' -mask_repeat RepeatMask' 
              . ' #expr( #password# ? " -dbpass " . #password# : "" )expr# ',
          },
      },
      {   -logic_name => 'export_transcript_sequences',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
              'cmd' => 'dump_genes.pl '
                . ' -dbuser #username# '
                . ' -dbname #dbname# '
                . ' -dbhost #host# '
                . ' -dbport #port# '
                . ' #expr( #password# ? " -dbpass " . #password# : "" )expr# '
                . ' -file #gene_sequences_file# '
                . ' -cdna -stable_id',
          },
          -rc_name   => '4Gb_job',
          -flow_into => {
              MEMLIMIT => 'export_transcript_sequences_himem',
          },
      },
      {   -logic_name => 'export_transcript_sequences_himem',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -rc_name   => '16Gb_job',
          -parameters => {
              'cmd' => 'dump_genes.pl '
                . ' -dbuser #username# '
                . ' -dbname #dbname# '
                . ' -dbhost #host# '
                . ' -dbport #port# '
                . ' #expr( #password# ? " -dbpass " . #password# : "" )expr# '
                . ' -file #gene_sequences_file# '
                . ' -cdna -stable_id',
          },
      },
      {
        -logic_name  => 'export_unmapped_sequences',
        -module      => 'Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::DumpUnmappedSeqs',
      },
    ];
}

1;
