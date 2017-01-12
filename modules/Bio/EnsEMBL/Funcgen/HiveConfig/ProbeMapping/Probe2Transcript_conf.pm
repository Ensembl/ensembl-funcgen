package Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::Probe2Transcript_conf;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Funcgen::HiveConfig::ProbeMapping::Base');

sub pipeline_analyses {
    my $self = shift;
    
    my $transcript_utr_file                            = '#tempdir#/unannotated_utrs.pl';
    my $extended_transcripts_file                      = '#tempdir#/extended_transcripts.bed';
    my $sorted_extended_transcripts_file               = '#tempdir#/extended_transcripts.sorted.bed';
    my $flanks_file                                    = '#tempdir#/flanks.pl';
    my $ungrouped_probe_feature_file                   = '#tempdir#/ungrouped_probe_features.bed';
    my $sorted_probe_feature_file                      = '#tempdir#/sorted_probe_features.bed';
    my $transcript_probe_features_overlaps_file        = '#tempdir#/transcript_probe_features_overlaps.bed';
    my $sorted_transcript_probe_features_overlaps_file = '#tempdir#/transcript_probe_features_overlaps.sorted.bed';
    my $arrays_per_object_file                         = '#tempdir#/arrays_per_object.pl';
    my $probeset_sizes_file                            = '#tempdir#/probeset_sizes.pl';
    my $object_names_file                              = '#tempdir#/object_names.pl';
    my $transcript_info_file                           = '#tempdir#/transcript_info.pl';
    my $probeset_transcript_assignments                = '#tempdir#/probeset_transcript_assignments.pl';
    my $probeset_transcript_rejections                 = '#tempdir#/probeset_transcript_rejections.pl';
    my $probe_transcript_assignments                   = '#tempdir#/probe_transcript_assignments.pl';

    return [
      {
          -logic_name  => 'start_probe2transcript',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -flow_into => {
              'MAIN->A' => [ 'calculate_utrs', 'export_probe_features_to_bed', 'calculate_arrays_per_object' ],
              'A->MAIN' => [ 'compute_transcript_probe_feature_overlaps'      ]
          },
      },
      {   -logic_name  => 'calculate_arrays_per_object',
          -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters  => {
              cmd => 
                  ' calculate_arrays_per_object.pl'
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
              cmd => "sort -u -k1,1 -k2,2n -k3,3n $extended_transcripts_file > $sorted_extended_transcripts_file"
          },
          -rc_name     => '8Gb_job',
      },
      
      mysql -N $(mysql-ens-reg-prod-1-ensadmin details mysql) homo_sapiens_core_86_38 -e 'select seq_region.name, seq_region.length from seq_region join seq_region_attrib using (seq_region_id) join attrib_type using (attrib_type_id) where code="toplevel" group by seq_region.name, seq_region.length order by seq_region.name' | sort -k1,1 > deleteme.txt
      
      
      {   -logic_name  => 'export_probe_features_to_bed',
          -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters  => {
              cmd => 
                'create_bedtools_genome_file.pl '
                . ' --registry #reg_conf#'
                . ' --species  #species#'
                . ' --file      ' . $bedtools_genome_file
          },
          -flow_into => {
              MAIN => 'sort_probe_feature_file',
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
                . ' -u' <--- Problem!?!
                . ' -T #tempdir#/#species#'
                . ' --buffer-size=31G'
                . ' -k1,1'
                . ' -k2,2n'
                . ' -k3,3n'
                .  ' --output=' . $sorted_probe_feature_file
                . ' ' . $ungrouped_probe_feature_file
          },
          -rc_name     => 'parallel_sort',
      },
      {   -logic_name  => 'compute_transcript_probe_feature_overlaps',
          -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters  => {
              cmd => 
                'bedtools'
                . ' intersect'
                . ' -sorted'
                . ' -wa -wb -a '.$extended_transcripts_file.' -b ' . $sorted_probe_feature_file
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
                . ' -u'
                . ' -T #tempdir#'
                . ' --buffer-size=31G'
                . ' -k4,4'
                . ' ' . $transcript_probe_features_overlaps_file
                .  ' > '
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
          },
          -flow_into => {
              MAIN => 'compute_hits',
          },
          -rc_name     => '4Gb_job',
      },
      {   -logic_name  => 'compute_hits',
          -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
          -parameters  => {
              cmd => 
                  'compute_hits.pl'
                . ' --probeset_sizes_file             ' . $probeset_sizes_file
                . ' --transcript_info_file            ' . $transcript_info_file
                . ' --probeset_transcript_assignments ' . $probeset_transcript_assignments
                . ' --probeset_transcript_rejections  ' . $probeset_transcript_rejections
                . ' --probe_transcript_assignments    ' . $probe_transcript_assignments
          },
      },


    ]
}

sub resource_classes {
    my ($self) = @_;
    return {
        %{$self->SUPER::resource_classes},
        'parallel_sort'     => {'LSF' => '-q production-rh7 -n 16 -M32000 -R"select[mem>32000] rusage[mem=32000]"' },
    };
}

1;
