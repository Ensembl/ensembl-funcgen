package Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::RunProbeToTranscriptHealthchecks_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::Base');

sub pipeline_analyses {
    my $self = shift;
    return [
      {
          -logic_name  => 'start_probe_to_transcript_healthchecks',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -flow_into => {
              MAIN => 'hc_methylation_arrays_not_mapped_to_transcripts'
          },
      },
      {
          -logic_name  => 'hc_methylation_arrays_not_mapped_to_transcripts',
          -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SqlHealthcheck',
          -parameters => {
            db_conn       => 'funcgen:#species#',
            description   => 'Assert probes from methylation arrays have not been mapped to transcripts',
            query         => "
              select 
                probe_id 
              from 
                array 
                join array_chip using (array_id) 
                join probe using (array_chip_id) 
                join probe_transcript using (probe_id) 
              where 
                format='METHYLATION' limit 1
            ",
            expected_size => '0'
          },
       -flow_into => {
           MAIN => 'hc_no_promiscuous_probes_mapped',
       },
      },
      {
          -logic_name  => 'hc_no_promiscuous_probes_mapped',
          -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SqlHealthcheck',
          -parameters => {
            db_conn       => 'funcgen:#species#',
            description   => 'Assert probes that are considered promiscuous have not been mapped to transcripts',
            query         => "
              select
                distinct probe_transcript.probe_id
              from
                unmapped_object
                join unmapped_reason on (
                  unmapped_object.unmapped_reason_id=unmapped_reason.unmapped_reason_id 
                  and summary_description = 'Promiscuous probe'
                )
                join probe_transcript on (probe_id = ensembl_id and ensembl_object_type = 'Probe')
              limit 1
            ",
            expected_size => '0'
          },
       -flow_into => {
           MAIN => 'hc_no_promiscuous_probesets_mapped',
       },
      },
      {
          -logic_name  => 'hc_no_promiscuous_probesets_mapped',
          -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SqlHealthcheck',
          -parameters => {
            db_conn       => 'funcgen:#species#',
            description   => 'Assert probe sets that are considered promiscuous have not been mapped to transcripts',
            query         => "
              select
                distinct probe_set_transcript.probe_set_id
              from
                unmapped_object
                join unmapped_reason on (
                  unmapped_object.unmapped_reason_id=unmapped_reason.unmapped_reason_id 
                  and summary_description = 'Promiscuous ProbeSet'
                )
                join probe_set_transcript on (probe_set_id = ensembl_id and ensembl_object_type = 'ProbeSet')
              limit 1
            ",
            expected_size => '0'
          },
          -flow_into => {
              MAIN => 'hc_probe_features_from_transcripts_accounted_for'
          },
      },
      {
          -logic_name  => 'hc_probe_features_from_transcripts_accounted_for',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -flow_into => {
              'MAIN->A' => 'hc_probe_features_from_transcripts_accounted_for_prepare',
              'A->MAIN' => 'hc_probe_features_from_transcripts_accounted_for_cleanup',
          },
      },
      {
          -logic_name  => 'hc_probe_features_from_transcripts_accounted_for_prepare',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
          -parameters => {
              sql     => [

                # Make sure the temporary tables don't already exist from a previous run.
                #
                'drop table if exists temp_probe_ids_with_transcript_matches',
                'drop table if exists temp_mapped_to_transcript',
                'drop table if exists temp_unmapped_by_probe_feature',
                'drop table if exists temp_probes_accounted_for',
                'drop table if exists temp_problem_probes',
                'drop table if exists temp_unmapped_promiscuous',

                # Collect all probe ids and the transcript on which they have probe features.
                # 
                # These probes should be either 
                #   - mapped to a transcript or 
                #   - the probe was classified as promiscuous or
                #   - the probe has at least one probe feature that was 
                #     rejected from being mapped to the transcript.
                #
                '
                  create table temp_probe_ids_with_transcript_matches as 
                    select 
                      distinct probe_id,
                      hit_id as  stable_id
                    from probe_feature where source="transcript";
                ',

                # Collect the probes that have been mapped to a transcript 
                # and the transcript they have been mapped to.
                #
                '
                  create table temp_mapped_to_transcript as 
                    select 
                      distinct probe_id, stable_id
                    from 
                      probe_transcript;
                ',

                # Collect the probes that have probe features that were 
                # rejected together with the transcript with which a mapping
                # was rejected.
                #
                '
                  create table temp_unmapped_by_probe_feature as 
                    select 
                      distinct probe_id, identifier as stable_id 
                    from 
                      unmapped_object 
                      join probe_feature on (ensembl_id=probe_feature_id and ensembl_object_type="ProbeFeature")
                ',

                # Probes that were classified as promiscuous
                '
                  create table temp_unmapped_promiscuous as
                    select 
                      ensembl_id as probe_id 
                    from 
                      unmapped_object 
                    where 
                      ensembl_object_type="Probe"
                ',
                'create index temp_unmapped_promiscuous_idx on temp_unmapped_promiscuous(probe_id);',

                # All probes in temp_probe_ids_with_transcript_matches should 
                # be in one of the three tables:
                #
                # - temp_probe_ids_with_transcript_matches meaning it has been mapped 
                #
                # or one of
                #
                # - temp_unmapped_by_probe_feature 
                # - temp_unmapped_promiscuous
                #
                # if it hasn't. If the probe is present in one of these three 
                # tables, it is accounted for.
                #
                '
                  create table temp_probes_accounted_for as 
                    select * from temp_mapped_to_transcript 
                    union select * from temp_unmapped_by_probe_feature
                ',
                'create index temp_probes_accounted_for_idx on temp_probes_accounted_for(probe_id, stable_id);',

                '
                  create table temp_problem_probes as 
                    select 
                      temp_probe_ids_with_transcript_matches.*
                    from 
                      temp_probe_ids_with_transcript_matches
                      left join temp_probes_accounted_for using (probe_id, stable_id)
                      left join temp_unmapped_promiscuous on (temp_probe_ids_with_transcript_matches.probe_id=temp_unmapped_promiscuous.probe_id)
                    where
                      temp_probes_accounted_for.probe_id is null
                      and temp_unmapped_promiscuous.probe_id is null
                    ;
                 ',

              ],
              db_conn => 'funcgen:#species#',
          },
          -flow_into => {
              MAIN => 'hc_probe_features_from_transcripts_accounted_for_run',
          },
      },
      {
          -logic_name  => 'hc_probe_features_from_transcripts_accounted_for_run',
          -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SqlHealthcheck',
          -parameters => {
            db_conn       => 'funcgen:#species#',
            description   => 'Find the number of probe features from transcripts whose probes have neither been mapped to a transcript nor are stored of have probe features that are stored as unmapped objects.',
            query         => 'select * from temp_problem_probes',
            expected_size => '0'
          },
      },
      {
          -logic_name  => 'hc_probe_features_from_transcripts_accounted_for_cleanup',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
          -parameters => {
              sql     => [
                'drop table temp_probe_ids_with_transcript_matches',
                'drop table temp_mapped_to_transcript',
                'drop table temp_unmapped_by_probe_feature',
                'drop table temp_probes_accounted_for',
                'drop table temp_problem_probes',
                'drop table temp_unmapped_promiscuous'
              ],
              db_conn => 'funcgen:#species#',
          },
      },
    ];
}

1;
