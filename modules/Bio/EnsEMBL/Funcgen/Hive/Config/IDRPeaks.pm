=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 NAME

    Bio::EnsEMBL::Funcgen::Hive::Config::IDRPeaks;

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 CONTACT

    Please contact http://lists.ensembl.org/mailman/listinfo/dev mailing list with questions/suggestions.

=cut
package Bio::EnsEMBL::Funcgen::Hive::Config::IDRPeaks;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Funcgen::Hive::Config::BaseSequenceAnalysis');

sub pipeline_wide_parameters {
  my $self = shift;
  return {
    %{$self->SUPER::pipeline_wide_parameters},
    
    can_run_SWEmbl_R0005_replicate => 1,#'IDRPeaks',
    can_PreprocessIDR              => 1,#'IDRPeaks',
    can_DefineMergedDataSet        => 0, 
   };
}

sub pipeline_analyses {
  my $self = shift;

  return [
   @{$self->SUPER::pipeline_analyses}, #To pick up BaseSequenceAnalysis-DefineMergedOutputSet
    {
     -logic_name    => 'PermissiveSWEmbl',
     -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunPeaks',
     -rc_name       => 'normal_5GB_2cpu_monitored',
    },
    {
     -logic_name    => 'PreprocessIDR',
#      -module        => 'Bio::EnsEMBL::Funcgen::Hive::PreprocessIDR',
#      -rc_name       => 'default',
#      -batch_size    => 30,
#      -parameters    => { 
# 	permissive_peaks => $self->o('permissive_peaks') 
#       },
     -flow_into => {
       '2->A' => 'RunIDR',
       'A->3' => 'PostProcessIDRReplicates', 
      }, 
    },
    {
      -logic_name    => 'RunIDR',
      -module        => 'Bio::EnsEMBL::Funcgen::Hive::RunIDR',
      -rc_name    => 'normal_2GB',
      -flow_into => {
	  2 => {
	    ':////accu?idr_peak_counts=[]' => {
	      'idr_peak_counts' => '#idr_peak_counts#'
	    },
	  },
# 	  2 => '?accu_name=idr_peak_counts&accu_address=[#idr_peak_counts#]'
      }
    },
    {
     -logic_name    => 'PostProcessIDRReplicates',
     -module        => 'Bio::EnsEMBL::Funcgen::Hive::PostprocessIDR',
     -rc_name    => 'default',
     -batch_size => 10,
     -flow_into => {
       2 => 'MergeIDRReplicateAlignments',
      }, 
    },
    {
      -logic_name    => 'MergeIDRReplicateAlignments',
      -module        => 'Bio::EnsEMBL::Funcgen::Hive::MergeIDRReplicateAlignments',
      -rc_name => 'default',
      -flow_into => { 
	2 => 'FixReplicateResultSetsExperimentIds'
      },
    },
    {
     -logic_name => 'FixReplicateResultSetsExperimentIds',
     -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SqlCmd',
     -parameters => {
	#
	# Sets the experiment id for the current result set.
	#
	# This should be set when creating the result set, but it is not. 
	# Until this is fixed in the api we do it here in an extra step.
	#
	sql => qq(
	  update result_set, epigenome, experiment 
	  set result_set.experiment_id = experiment.experiment_id 
	  where epigenome.epigenome_id=experiment.epigenome_id and epigenome.epigenome_id=result_set.epigenome_id and experiment.feature_type_id=result_set.feature_type_id
	  and result_set_id = #dbID#
	  ),
	  db_conn => 'funcgen:#species#',
	},
#       -meadow     => 'LOCAL',
      -flow_into => {
	1 => 'DefineMergedDataSet',
      },
    },
  ];
}

1;
