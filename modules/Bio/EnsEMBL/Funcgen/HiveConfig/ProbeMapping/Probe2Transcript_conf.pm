package Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::Probe2Transcript_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::Base');

sub pipeline_analyses {
    my $self = shift;
    
    my $probe2transcript_temp_dir                = '#tempdir#/#species#/probe2transcript';
    my $probe2transcript_array_specific_temp_dir = $probe2transcript_temp_dir . '/#array_name#';

    my $transcript_utr_file                            = $probe2transcript_temp_dir . '/unannotated_utrs.pl';
    my $extended_transcripts_file                      = $probe2transcript_temp_dir . '/extended_transcripts.bed';
    my $sorted_extended_transcripts_file               = $probe2transcript_temp_dir . '/extended_transcripts.sorted.bed';
    my $flanks_file                                    = $probe2transcript_temp_dir . '/flanks.pl';
    my $bedtools_genome_file                           = $probe2transcript_temp_dir . '/bedtools_genome_file';
    my $ungrouped_probe_feature_file                   = $probe2transcript_temp_dir . '/ungrouped_probe_features.bed';
    my $sorted_probe_feature_file                      = $probe2transcript_temp_dir . '/sorted_probe_features.bed';
    my $transcript_probe_features_overlaps_file        = $probe2transcript_temp_dir . '/transcript_probe_features_overlaps.bed';
    my $sorted_transcript_probe_features_overlaps_file = $probe2transcript_temp_dir . '/transcript_probe_features_overlaps.sorted.bed';
    my $arrays_per_object_file                         = $probe2transcript_temp_dir . '/arrays_per_object.pl';
    my $probeset_sizes_file                            = $probe2transcript_temp_dir . '/probeset_sizes.pl';
    my $object_names_file                              = $probe2transcript_temp_dir . '/object_names.pl';
    my $transcript_info_file                           = $probe2transcript_temp_dir . '/transcript_info.pl';
    my $probe_feature_transcript_rejection_file        = $probe2transcript_temp_dir . '/rejected_probe_features.pl';

    my $probeset_transcript_hits_by_array_file         = $probe2transcript_array_specific_temp_dir . '/probeset_transcript_hits_file.pl';
    my $probe_transcript_hits_by_array_file            = $probe2transcript_array_specific_temp_dir . '/probe_transcript_hits_file.pl';
    my $probeset_to_transcript_by_array_file           = $probe2transcript_array_specific_temp_dir . '/probeset_to_transcript_file.pl';
    my $probe_transcript_assignments_by_array_file     = $probe2transcript_array_specific_temp_dir . '/probes_transcript_assignments.tsv';
    my $probeset_transcript_assignments_by_array       = $probe2transcript_array_specific_temp_dir . '/probeset_transcript_assignments.tsv';
    my $probeset_rejections_file                       = $probe2transcript_array_specific_temp_dir . '/rejected_probesets.pl';
    
    my $probe_feature_transcript_assignment_file       = $probe2transcript_temp_dir . '/probe_feature_transcript_assignments.tsv';
    
    my $probeset_transcript_rejections                 = $probe2transcript_temp_dir . '/probeset_transcript_rejections.tsv';
    my $probe_transcript_assignments                   = $probe2transcript_temp_dir . '/probe_transcript_assignments.tsv';

    return [
      {
          -logic_name  => 'start_probe2transcript',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -flow_into => {
              MAIN => 'truncate_p2t_tables'
          },
      },
      {
          -logic_name  => 'truncate_p2t_tables',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
          -parameters => {
              sql     => [
                "truncate probe_transcript;",
                "truncate probe_set_transcript;",
                "truncate probe_feature_transcript;",
                "truncate unmapped_object;",
                "truncate unmapped_reason;",
              ],
              db_conn => 'funcgen:#species#',
          },
          -flow_into => {
              MAIN => 'switch_p2t_tables_to_innodb',
          },
      },
      {
          -logic_name  => 'switch_p2t_tables_to_innodb',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
          -parameters => {
              sql     => [
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
              MAIN => 'p2t_mk_tempdir',
          },
      },
      {   -logic_name  => 'p2t_mk_tempdir',
          -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters  => {
              cmd => 'mkdir -p ' . $probe2transcript_temp_dir
          },
          -flow_into => {
              'MAIN->A' => [ 
                'calculate_utrs', 
                'export_probe_features_to_bed', 
                'create_bedtools_genome_file', 
                'calculate_arrays_per_object' 
              ],
              'A->MAIN' => 'compute_transcript_probe_feature_overlaps'
          },
      },
      {   -logic_name  => 'calculate_arrays_per_object',
          -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters  => {
              cmd => 
                  'calculate_arrays_per_object.pl'
                . ' --registry #reg_conf#'
                . ' --species  #species#'
                . ' --arrays_per_object_file ' . $arrays_per_object_file
                . ' --probeset_sizes_file    ' . $probeset_sizes_file
                . ' --object_names_file      ' . $object_names_file
          },
      },
      {   -logic_name  => 'calculate_utrs',
          -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters  => {
              cmd => 'calculate_utrs.pl'
                . ' --registry #reg_conf#'
                . ' --species  #species#'
                . ' --transcript_utr_file ' . $transcript_utr_file,
          },
          -flow_into => {
              MAIN => 'write_extended_transcripts_into_file',
          },
          -rc_name     => '8Gb_job',
      },
      {   -logic_name  => 'write_extended_transcripts_into_file',
          -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters  => {
              cmd => 'write_extended_transcripts_into_file.pl '
                . ' --registry #reg_conf#'
                . ' --species  #species#'
                . ' --unannotated_utrs '               . $transcript_utr_file
                . ' --flanks_outputfile '              . $flanks_file
                . ' --extended_transcript_outputfile ' . $extended_transcripts_file
          },
          -flow_into => {
              MAIN => 'sort_extended_transcripts_file',
          },
          -rc_name     => '8Gb_job',
      },
      {   -logic_name  => 'sort_extended_transcripts_file',
          -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters  => {
#               cmd => 'bedSort ' . $extended_transcripts_file . ' ' . $extended_transcripts_file
              # Can't use bedSort or the chromosome names won't be consistent with $sorted_probe_feature_file
              #
              cmd => "sort -k1,1 -k2,2n -k3,3n $extended_transcripts_file | uniq > $sorted_extended_transcripts_file"
          },
          -rc_name     => '8Gb_job',
      },
      {   -logic_name  => 'create_bedtools_genome_file',
          -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters  => {
              cmd => 
                'create_bedtools_genome_file.pl '
                . ' --registry #reg_conf#'
                . ' --species  #species#'
                . ' --bedtools_genome_file ' . $bedtools_genome_file
          },
      },
      {   -logic_name  => 'export_probe_features_to_bed',
          -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters  => {
              cmd => 
                'export_probe_features_to_bed.pl '
                . ' --registry #reg_conf#'
                . ' --species  #species#'
                . ' --file      ' . $ungrouped_probe_feature_file
          },
          -flow_into => {
              MAIN => 'sort_probe_feature_file',
          },
      },
      {   -logic_name  => 'sort_probe_feature_file',
          -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters  => {
              # Monitor on node using:
              # ps -C sort -L -o pcpu
              #
              cmd => 
                'sort'
                . ' --parallel=16'
                . ' -T ' . $probe2transcript_temp_dir
                . ' --buffer-size=31G'
                . ' -k1,1'
                . ' -k2,2n'
                . ' -k3,3n'
                . ' ' . $ungrouped_probe_feature_file
                # Sort has a -u option, but it doesn't seem to look at the 
                # entire line, only the columns from the -k parameters.
                #
                . ' | uniq '
                . ' > ' . $sorted_probe_feature_file
          },
          -rc_name     => 'parallel_sort',
      },
      {   -logic_name  => 'compute_transcript_probe_feature_overlaps',
          -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters  => {
              cmd => 
                'bedtools'
                . ' intersect'
                . ' -g ' . $bedtools_genome_file
                . ' -sorted'
                . ' -wa -wb -a '.$sorted_extended_transcripts_file .' -b ' . $sorted_probe_feature_file
                . ' > ' . $transcript_probe_features_overlaps_file
          },
          -flow_into => {
              MAIN => 'sort_overlaps_by_transcript_stable_id',
          },
      },
      {   -logic_name  => 'sort_overlaps_by_transcript_stable_id',
          -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters  => {
              cmd => 
                'sort'
                . ' --parallel=16'
                . ' -T ' . $probe2transcript_temp_dir
                . ' --buffer-size=31G'
                . ' -k4,4'
                . ' ' . $transcript_probe_features_overlaps_file
                .  ' | uniq > '
                . $sorted_transcript_probe_features_overlaps_file
          },
          -rc_name     => 'parallel_sort',
          -flow_into => {
              MAIN => 'examine_transcript',
          },
      },
      {   -logic_name  => 'examine_transcript',
          -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters  => {
              cmd => 
                  'examine_transcript.pl'
                . ' --registry #reg_conf#'
                . ' --species  #species#'
                . ' --transcript_probe_features_overlaps ' . $sorted_transcript_probe_features_overlaps_file
                . ' --flanks_file                        ' . $flanks_file
                . ' --transcript_utr_file                ' . $transcript_utr_file 
                . ' --transcript_info_file               ' . $transcript_info_file
                . ' --probe_feature_transcript_assignments_file ' . $probe_feature_transcript_assignment_file
                . ' --probe_feature_transcript_rejection_file '   . $probe_feature_transcript_rejection_file
          },
          -flow_into => {
              MAIN => [
                'compute_probeset_transcript_assignments_per_array',
                'load_probe_feature_to_transcript_assignments',
                'load_probe_feature_to_transcript_rejections',
              ]
          },
          -rc_name     => '4Gb_job',
      },
        {   -logic_name  => 'load_probe_feature_to_transcript_rejections',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 
                    'load_probeset_to_transcript_rejections.pl '
                  . ' --registry    #reg_conf#'
                  . ' --species     #species#'
                  . ' --analysis_logic_name  probe2transcript'
                  . ' --probeset_rejections_file ' . $probe_feature_transcript_rejection_file
            },
        },

        {   -logic_name  => 'load_probe_feature_to_transcript_assignments',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 
                    'load_probe_feature_to_transcript_assignments.pl'
                  . ' --registry    #reg_conf#'
                  . ' --species     #species#'
                  . ' --array_name  foo'
                  . ' --probe_feature_transcript_assignments_file ' . $probe_feature_transcript_assignment_file
            },
        },
        {   -logic_name  => 'compute_probeset_transcript_assignments_per_array',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -parameters => {
                db_conn    => 'funcgen:#species#',
                inputquery => 'select "#species#" as species, array.name as array_name from array where format!="METHYLATION"',
            },
            -flow_into => {
               2 => 'p2t_mk_arrays_tempdir',
            },
        },
        {   -logic_name  => 'p2t_mk_arrays_tempdir',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 'mkdir -p ' . $probe2transcript_array_specific_temp_dir
            },
            -flow_into => {
                MAIN => [
                  'compute_probeset_transcript_hits',
                  'compute_probe_transcript_hits'
                ]
            },
        },
        {   -logic_name  => 'compute_probeset_transcript_hits',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 
                    'compute_probeset_transcript_hits.pl'
                  . '  --array_name #array_name#'
                  . '  --transcript_info_file '          . $transcript_info_file
                  . '  --probeset_transcript_hits_file ' . $probeset_transcript_hits_by_array_file
            },
            -flow_into => {
                MAIN     => 'create_probeset_to_transcript_descriptions',
                MEMLIMIT => 'compute_probeset_transcript_hits_32gb',
            },
        },
        {   -logic_name  => 'compute_probeset_transcript_hits_32gb',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 
                    'compute_probeset_transcript_hits.pl'
                  . '  --array_name #array_name#'
                  . '  --transcript_info_file       ' . $transcript_info_file
                  . '  --probeset_transcript_hits_file ' . $probeset_transcript_hits_by_array_file
            },
            -flow_into => {
                MAIN => 'create_probeset_to_transcript_descriptions',
                MEMLIMIT => 'compute_probeset_transcript_hits_64gb',
            },
            -rc_name     => '32Gb_job',
        },
        {   -logic_name  => 'compute_probeset_transcript_hits_64gb',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 
                    'compute_probeset_transcript_hits.pl'
                  . '  --array_name #array_name#'
                  . '  --transcript_info_file       ' . $transcript_info_file
                  . '  --probeset_transcript_hits_file ' . $probeset_transcript_hits_by_array_file
            },
            -flow_into => {
                MAIN => 'create_probeset_to_transcript_descriptions',
            },
            -rc_name     => '64Gb_job',
        },
        {   -logic_name  => 'create_probeset_to_transcript_descriptions',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 
                    'create_probeset_to_transcript_descriptions.pl '
                  . '  --array_name #array_name#'
                  . '  --probeset_sizes                      '    . $probeset_sizes_file
                  . '  --probeset_transcript_hits_by_array_file ' . $probeset_transcript_hits_by_array_file
                  . '  --probeset_to_transcript_file         '    . $probeset_transcript_assignments_by_array
                  . '  --rejected_probesets_file ' . $probeset_rejections_file
            },
            -flow_into => {
                MAIN => [
                  'load_probeset_to_transcript_assignments',
                  'load_probeset_to_transcript_rejections'
                ]
            },
        },
        {   -logic_name  => 'load_probeset_to_transcript_rejections',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 
                    'load_probeset_to_transcript_rejections.pl '
                  . ' --registry    #reg_conf#'
                  . ' --species     #species#'
                  . ' --analysis_logic_name  ProbeAlign_transcript'
                  . ' --probeset_rejections_file ' . $probeset_rejections_file
            },
        },
        {   -logic_name  => 'load_probeset_to_transcript_assignments',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 
                    'load_probeset_to_transcript_assignments.pl'
                  . ' --registry    #reg_conf#'
                  . ' --species     #species#'
                  . ' --array_name  #array_name#'
                  . ' --probeset_transcript_assignments_file ' . $probeset_transcript_assignments_by_array
            },
        },
        {   -logic_name  => 'compute_probe_transcript_hits',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 
                    'compute_probe_transcript_hits.pl'
                  . '  --array_name #array_name#'
                  . '  --transcript_info_file       ' . $transcript_info_file
                  . '  --probe_transcript_hits_file ' . $probe_transcript_hits_by_array_file
            },
            -flow_into => {
                MAIN     => 'create_probe_to_transcript_descriptions',
                MEMLIMIT => 'compute_probe_transcript_hits_32gb',
            },
        },
        {   -logic_name  => 'compute_probe_transcript_hits_32gb',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 
                    'compute_probe_transcript_hits.pl'
                  . '  --array_name #array_name#'
                  . '  --transcript_info_file       ' . $transcript_info_file
                  . '  --probe_transcript_hits_file ' . $probe_transcript_hits_by_array_file
            },
            -flow_into => {
                MAIN     => 'create_probe_to_transcript_descriptions',
                MEMLIMIT => 'compute_probe_transcript_hits_64gb',
            },
          -rc_name     => '32Gb_job',
        },
        {   -logic_name  => 'compute_probe_transcript_hits_64gb',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 
                    'compute_probe_transcript_hits.pl'
                  . '  --array_name #array_name#'
                  . '  --transcript_info_file       ' . $transcript_info_file
                  . '  --probe_transcript_hits_file ' . $probe_transcript_hits_by_array_file
            },
            -flow_into => {
                MAIN => 'create_probe_to_transcript_descriptions',
            },
          -rc_name     => '64Gb_job',
        },
        {   -logic_name  => 'create_probe_to_transcript_descriptions',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 
                    'create_probe_to_transcript_descriptions.pl'
                  . '  --array_name #array_name#'
                  . '  --probe_transcript_hits_by_array_file ' . $probe_transcript_hits_by_array_file
                  . '  --probe_to_transcript_file            ' . $probe_transcript_assignments_by_array_file
            },
            -flow_into => {
                MAIN => 'load_probe_to_transcript_assignments',
            },
        },
        {   -logic_name  => 'load_probe_to_transcript_assignments',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                cmd => 
                    'load_probe_to_transcript_assignments.pl'
                  . ' --registry    #reg_conf#'
                  . ' --species     #species#'
                  . ' --array_name  #array_name#'
                  . ' --probe_transcript_assignments_file ' . $probe_transcript_assignments_by_array_file
            },
        },
    ]
}

sub resource_classes {
    my ($self) = @_;
    return {
        %{$self->SUPER::resource_classes},
        'parallel_sort'     => { 'LSF' => '-q production-rh7 -n 16 -M32000 -R"select[mem>32000] rusage[mem=32000]"' },
    };
}

1;
