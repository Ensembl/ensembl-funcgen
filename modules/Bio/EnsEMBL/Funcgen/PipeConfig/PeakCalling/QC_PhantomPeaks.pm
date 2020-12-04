=pod 

=head1 NAME

    Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::QC_PhantomPeaks

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
package Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::QC_PhantomPeaks;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
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
              2 => 'qc_phantom_peaks_start',
            },
        },
        {   -logic_name => 'qc_phantom_peaks_start',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
            -flow_into => { 
                'MAIN->A' => 'qc_phantom_peaks_job_factory',
                'A->MAIN' => 'qc_phantom_peaks_done',
            },
        },
        {   -logic_name => 'qc_phantom_peaks_job_factory',
            -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::QcPhantomPeaksJobFactory',
            -parameters  => {
              tempdir => '#tempdir_peak_calling#/#species#/phantom_peaks',
            },
            -flow_into => { 
                2 => 'qc_run_phantom_peaks',
            },
        },
        {   -logic_name  => 'qc_run_phantom_peaks',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                  use_bash_pipefail => 1,
                  return_codes_2_branches => {
                    1 => 0
                  },
                  cmd => 
                      # Rscript does not search the path, so we use "which" to 
                      # do that. Also using single quotes to avoid interpolation
                      # of the dollar sign.
                      #
                      q( Rscript $(which run_spp_nodups.R)   )
                      # Overwrite plotfile, if one already exists
                      . qq(    -rf                             )
                      . qq(    -c=#bam_file#                   )
                      . qq(    -savp                           )
                      . qq(    -out=#phantom_peak_out_file#    )
                      . qq(    -odir=#phantom_peak_tempdir#    )
                      . qq(    -tmpdir=#phantom_peak_tempdir#  )
                      # Avoid the "Error: ignoring SIGPIPE signal error"
                      # from the EBI cluster
                      . qq(    &&                              )
                      # In case the job gets terminated for memlimit, this 
                      # ensures that the worker also dies. (or so we hope)
                      . qq(    sleep 30                        )
            },
          -rc_name    => '10Gb_job_2cpus',
          -flow_into  => { 
              MAIN     => 'qc_load_phantom_peaks',
              0        => 'qc_run_phantom_peaks_himem',
          },
        },
        {   -logic_name  => 'qc_run_phantom_peaks_himem',
            -module      => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters  => {
                  use_bash_pipefail => 1,
                  return_codes_2_branches => {
                    1 => 0
                  },
                  cmd => 
                  # Rscript does not search the path, so we use "which" to 
                  # do that. Also using single quotes to avoid interpolation
                  # of the dollar sign.
                  #
                  q( Rscript $(which run_spp_nodups.R)   )
                  # Overwrite plotfile, if one already exists
                  . qq(    -rf                             )
                  . qq(    -c=#bam_file#                   )
                  . qq(    -savp                           )
                  . qq(    -out=#phantom_peak_out_file#    )
                  . qq(    -odir=#phantom_peak_tempdir#    )
                  . qq(    -tmpdir=#phantom_peak_tempdir#  )
                  # Avoid the "Error: ignoring SIGPIPE signal error"
                  # from the EBI cluster
                  . qq(    &&                              )
                  # In case the job gets terminated for memlimit, this 
                  # ensures that the worker also dies. (or so we hope)
                  . qq(    sleep 30                        )
            },
	    -rc_name    => '32Gb_job_2cpus',
            -flow_into  => { 
              MAIN => 'qc_load_phantom_peaks',
              0    => 'qc_load_failed_phantom_peaks',
            },
        },
        {   -logic_name => 'qc_load_phantom_peaks',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd =>
                    qq( load_phantom_peak_file.pl                   )
                  . qq(    --alignment_name #alignment_name#        )
                  . qq(    --result_file    #phantom_peak_out_file# )
                  . qq(    --species        #species#               )
                  . qq(    --registry       #reg_conf#              )
                  . qq(    --failed         0                )
            },
        },
        {   -logic_name => 'qc_load_failed_phantom_peaks',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -parameters => {
                cmd =>
                    qq( load_phantom_peak_file.pl                   )
                  . qq(    --alignment_name #alignment_name#        )
                  . qq(    --result_file    #phantom_peak_out_file# )
                  . qq(    --species        #species#               )
                  . qq(    --registry       #reg_conf#              )
                  . qq(    --failed         1                )
            },
        },
        {   -logic_name => 'qc_phantom_peaks_done',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
        },
    ];
}

1;
