=pod 

=head1 NAME

    Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::QC_Chance

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 LICENSE

    Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
    Copyright [2016-2020] EMBL-European Bioinformatics Institute

    Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

         http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software distributed under the License
    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and limitations under the License.

=head1 CONTACT

=cut

package Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::QC_Chance;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
#use base ('Bio::EnsEMBL::Funcgen::Hive::Config::Base');
use base 'Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::Base';

sub pipeline_analyses {
    my ($self) = @_;
    return [
        {   -logic_name => 'start_alignment_qc',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into => { 
              1 => 'alignment_qc_seed_all_plans',
            },
        },
        {   -logic_name => 'alignment_qc_seed_all_plans',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::SeedAllExecutionPlans',
            -flow_into => { 
              2 => 'alignment_qc_seed_all_alignment',
            },
        },
        {   -logic_name => 'alignment_qc_seed_all_alignment',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::SeedAllAlignments',
            -flow_into => { 
              2 => 'qc_chance_start',
            },
        },
        {   -logic_name => 'qc_chance_start',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into => { 
              'MAIN->A' => 'QcChanceJobFactory',
              'A->MAIN' => 'qc_chance_done',
            },
        },
        {
          -logic_name    => 'QcChanceJobFactory',
          -module        => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::QcChanceJobFactory',
          -flow_into => { 
              2 => [
                  WHEN('not(#alignment_has_duplicates#)' => ['MkTempDir']),
              ]
          },
        },
        {   -logic_name => 'MkTempDir',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd => qq!mkdir -p #chance_tempdir#!,
            },
            -flow_into => { 
              'MAIN->A' => [ 'JobFactoryArgenrich', 'CreateChanceBins'],
              'A->MAIN' => 'RunArgenrich',
            },
        },
        {   -logic_name => 'CreateChanceBins',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd => qq(create_chance_bins.pl )
                . qq(  --species                 #species#                     )
                . qq(  --epigenome_gender        #epigenome_gender#            )
                . qq(  --assembly                #assembly#                    )
                . qq(  --outputfile              #chance_tempdir#/#chance_bin_file#   )
                . qq(  --reference_data_root_dir #reference_data_root_dir#     )
            },
        },
        {   -logic_name => 'JobFactoryArgenrich',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -flow_into => {
                2 => 'CpToTemp',
            },
        },
        {   -logic_name => 'CpToTemp',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => { 
                  cmd => qq!rm -f #chance_tempdir#/#file# ; ln -s #sourcedir#/#file# #chance_tempdir#!,
            },
            -flow_into => { MAIN => 'IndexBam' },
        },
        {   -logic_name => 'IndexBam',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => { 
                cmd => qq!samtools index #chance_tempdir#/#file#!,
            },
            -analysis_capacity => 50,
            -flow_into => { 
              MAIN     => 'CountReads',
              MEMLIMIT => 'IndexBam_himem',
            },
        },
        {   -logic_name => 'IndexBam_himem',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => { 
                cmd => qq!samtools index #chance_tempdir#/#file#!,
            },
            -rc_name   => '8Gb_job',
            -flow_into => { MAIN => 'CountReads' },
        },
        {   -logic_name => 'CountReads',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
            -parameters => { 
                inputcmd        => "samtools view -c #chance_tempdir#/#file#",
                column_names    => [ 'read_count' ],
            },
            -flow_into => {
                2 => [
                    '?accu_name=read_count&accu_address={kind}',
                    '?accu_name=file&accu_address={kind}',
                ],
             },
        },
        {   -logic_name => 'RunArgenrich',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -max_retry_count => 1,
            -parameters => {
                cmd => qq(timeout --kill-after=0 12h argenrich_with_labels_and_rerunnable.R --args plot=TRUE outdir=#chance_tempdir# )
                  . qq(    ipsz=#expr( #read_count#->{"signal"}            )expr# )
                  . qq( inputsz=#expr( #read_count#->{"control"}           )expr# )
                  . qq(      ip=#chance_tempdir#/#expr( #file#->{"signal"}        )expr# )
                  . qq(   input=#chance_tempdir#/#expr( #file#->{"control"}       )expr# )
                  . qq(    bins=#chance_tempdir#/#chance_bin_file#                       )

                  # This ends up in #chance_tempdir#, because of the parameter 
                  # "outdir=#chance_tempdir#" set further above.
                  #
                  . qq( outfile=#argenrich_outfile#),
                  return_codes_2_branches => {
                    1   => 2,
                    
                    # "Unable to iterate to region within BAM."
                    124 => 2
                  },
            },
            -flow_into => {
                MAIN => 'load_chance',
                2    => 'load_chance_failed', 
            },
            -rc_name   => '16Gb_job',
        },
        {   -logic_name => 'load_chance',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd => qq(load_argenrich_qc_file.pl   )
                . qq( --argenrich_file        #chance_tempdir#/#argenrich_outfile#     )
                . qq( --signal   #signal_alignment#       )
                . qq( --control  #control_alignment#      )
                . qq( --user     #tracking_db_user#       )
                . qq( --pass     #tracking_db_pass#       )
                . qq( --port     #tracking_db_port#       )
                . qq( --host     #tracking_db_host#       )
                . qq( --dbname   #tracking_db_name#       )
                . qq( --work_dir #chance_tempdir#                )
                . qq( --species  #species#                )
                . qq( --failed   0                 )
            },
        },
        {   -logic_name => 'load_chance_failed',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd => qq(load_argenrich_qc_file.pl   )
                . qq( --argenrich_file        #chance_tempdir#/#argenrich_outfile#     )
                . qq( --signal   #signal_alignment#       )
                . qq( --control  #control_alignment#      )
                . qq( --user     #tracking_db_user#       )
                . qq( --pass     #tracking_db_pass#       )
                . qq( --port     #tracking_db_port#       )
                . qq( --host     #tracking_db_host#       )
                . qq( --dbname   #tracking_db_name#       )
                . qq( --work_dir #chance_tempdir#                )
                . qq( --species  #species#                )
                . qq( --failed   1                 )
            },
        },
        {   -logic_name => 'qc_chance_done',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
    ];
}

1;
