package Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::RunAlignHealthchecks_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::Base');

my $max_allowed_hits_per_probe = 100;

sub pipeline_analyses {
    my $self = shift;
    
    my $max_probe_features_to_check = 10000;
    
    return [
      {
          -logic_name  => 'start_align_healthchecks',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -flow_into => {
              MAIN => 'check_duplicate_probe_features'
          },
      },
      {   -logic_name        => 'check_duplicate_probe_features',
          -module            => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
              cmd => '
                delete_redundant_probe_features.pl \
                  --registry #reg_conf# \
                  --species #species# \
                  --only_test 1
              ',
          },
          -flow_into => {
              MAIN => 'check_gapped_probe_features_from_transcript_matches'
          },
      },
      {
          -logic_name  => 'check_gapped_probe_features_from_transcript_matches',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters => {
              cmd => 'check_probe_feature_sequences.pl '
                . ' --registry #reg_conf#'
                . ' --species  #species#'
                . ' --logic_name ProbeAlign_transcript '
                . ' --check_probe_features_with_nontrivial_cigar_lines 1'
                . ' --max_check ' . $max_probe_features_to_check
          },
          -flow_into => {
              MAIN => 'check_ungapped_probe_features_from_transcript_matches',
          },
        },
      {
          -logic_name  => 'check_ungapped_probe_features_from_transcript_matches',
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
          -flow_into => {
              MAIN => 'hc_no_probe_features_from_known_promiscuous_probes',
          },
      },
      {
          -logic_name  => 'hc_no_probe_features_from_known_promiscuous_probes',
          -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SqlHealthcheck',
          -parameters => {
            db_conn       => 'funcgen:#species#',
            description   => 'Assert that there are no probe features from probes that were classified as promiscuous',
            query         => "
              select
                distinct probe_feature.probe_feature_id
              from
                unmapped_object
                join unmapped_reason on (
                  unmapped_object.unmapped_reason_id=unmapped_reason.unmapped_reason_id and summary_description = 'Promiscuous probe'
                )
                join probe_feature on (probe_id = ensembl_id and ensembl_object_type = 'Probe')
              limit 1
            ",
            expected_size => '0'
          },
          -flow_into => {
              MAIN => 'hc_no_probe_features_from_promiscuous_probes',
          },
      },
      {
          -logic_name  => 'hc_no_probe_features_from_promiscuous_probes',
          -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SqlHealthcheck',
          -parameters => {
            db_conn       => 'funcgen:#species#',
            description   => 'Assert that there are no probe features from probes that were classified as promiscuous',
            query         => "select probe_id, count(probe_feature_id) c from probe_feature group by probe_id having c > $max_allowed_hits_per_probe order by c desc",
            expected_size => '0'
          },
      },
    ];
}

1;
