=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

    Bio::EnsEMBL::Funcgen::Hive::Config::DefineMergedDataSets;

=head1 SYNOPSIS


=head1 DESCRIPTION


=head1 CONTACT

    Please contact http://lists.ensembl.org/mailman/listinfo/dev mailing list with questions/suggestions.

=cut

package Bio::EnsEMBL::Funcgen::Hive::Config::DefineMergedDataSets;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Funcgen::Hive::Config::BaseSequenceAnalysis');

sub pipeline_wide_parameters {
  my $self = shift;
  return {
    %{$self->SUPER::pipeline_wide_parameters}, 
    can_DefineMergedDataSet => 1,
   };
}

sub pipeline_analyses {
  my $self = shift;

  return [
    {
      -logic_name => 'DefineMergedDataSet', 
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::DefineDataSet',
      -parameters => {
	default_feature_set_analyses => $self->o('default_peak_analyses'),
	feature_set_analysis_type    => 'peak',
      },
      -flow_into => {
	1 => [ 'PreprocessAlignments' ],
      },
      -analysis_capacity => 100,
      -rc_name           => 'default',
      #None of these shoudl take > 2-3 mins unless there is some rolling back to do
      #But this will only ever be deleting annotated_feature records
      #and maybe some states or updating sets
      #Keep this fairly low, so we a getting the parallel compute running asap.
      -batch_size        => 10, 
    },
    {
      -logic_name => 'PreprocessAlignments',
      # This is basically making sure the input file is sorted wrt genomic locations
      -module     => 'Bio::EnsEMBL::Funcgen::Hive::CollectionWriter',
      #Need to maintain this here as will not be updated by -analysis_topup
      -parameters => {
      feature_formats => ['bam', 'bed'],
      },
      -analysis_capacity => 100,
      -rc_name => 'normal_high_mem_2cpu',
    },
  ];
}

1;
