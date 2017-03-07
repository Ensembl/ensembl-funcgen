package Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::RunHealthchecks_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::Base');

sub pipeline_analyses {
    my $self = shift;
    
    my $max_probe_features_to_check = 10000;
    
    return [
      {
          -logic_name  => 'start_healthchecks',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -flow_into => {
              MAIN => 'check_probe_feature_sequences_from_transcripts_matches'
          },
      },
      {
          -logic_name  => 'check_probe_feature_sequences_from_transcripts_matches',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
              cmd => 'check_probe_feature_sequences.pl '
                . ' --registry #reg_conf#'
                . ' --species  #species#'
                . ' --logic_name ProbeAlign_transcript '
                . ' --max_check ' . $max_probe_features_to_check
          },
          -flow_into => {
              MAIN => 'check_probe_feature_sequences_from_genomic_matches',
          },
        },
      {
          -logic_name  => 'check_probe_feature_sequences_from_genomic_matches',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
              cmd => 'check_probe_feature_sequences.pl '
                . ' --registry #reg_conf#'
                . ' --species  #species#'
                . ' --logic_name ProbeAlign_genomic '
                . ' --max_check ' . $max_probe_features_to_check
          },
        },
    ];
}

1;
