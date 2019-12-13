=pod 

=head1 NAME

    Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::QC_ProportionOfReadsInPeaks

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
package Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::QC_ProportionOfReadsInPeaks;

use strict;
use warnings;

use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;
use base 'Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::Base';

sub pipeline_wide_parameters {
    my ($self) = @_;
    return {
      %{$self->SUPER::pipeline_wide_parameters},

      reg_conf => $self->o('reg_conf'),
    };
}

sub pipeline_analyses {
    my $self = shift;
    return [
      {
          -logic_name  => 'start_frip',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -flow_into => { 
            1 => 'qc_frip_seed_all_execution_plans',
          },
      },
      {
          -logic_name  => 'qc_frip_seed_all_execution_plans',
          -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::SeedAllExecutionPlans',
          -flow_into => {
            2 => 'seed_all_peak_callings_from_experiment',
          },
      },
      {
          -logic_name  => 'seed_all_peak_callings_from_experiment',
          -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::SeedAllPeakCallingsFromExperiment',
          -flow_into => { 
            2 => 'start_qc_proportion_of_reads_in_peaks',
          },
      },
      {   -logic_name => 'start_qc_proportion_of_reads_in_peaks',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
          -flow_into => { 
            'MAIN->A' => 'qc_prepare_proportion_of_reads_in_peaks',
            'A->MAIN' => 'done_qc_proportion_of_reads_in_peaks',
          },
      },
      {   -logic_name => 'qc_prepare_proportion_of_reads_in_peaks',
          -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::QcProportionOfReadsInPeaksJobFactory',
          -parameters => {
            tempdir => '#tempdir_peak_calling#/#species#/frip',
          },
          -flow_into => { 
            2 => 'qc_compute_proportion_of_reads_in_peaks',
          },
      },
      {   -logic_name => 'qc_compute_proportion_of_reads_in_peaks',
          -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::QcComputeProportionOfReadsInPeaks',
          -analysis_capacity => 50,
          -flow_into => { 
            2 => 'qc_store_proportion_of_reads_in_peaks',
            MEMLIMIT => 'qc_compute_proportion_of_reads_in_peaks_himem',
          },
      },
      {   -logic_name => 'qc_compute_proportion_of_reads_in_peaks_himem',
          -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::QcComputeProportionOfReadsInPeaks',
          -analysis_capacity => 50,
          -rc_name => '8Gb_job',
          -flow_into => { 
            2 => 'qc_store_proportion_of_reads_in_peaks',
          },
      },
      
      {   -logic_name => 'qc_store_proportion_of_reads_in_peaks',
          -module     => 'Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::QcStoreProportionOfReadsInPeaks',
          -analysis_capacity => 1,
      },
      {   -logic_name => 'done_qc_proportion_of_reads_in_peaks',
          -module     => 'Bio::EnsEMBL::Hive::RunnableDB::Dummy',
      },
    ];
}

1;
