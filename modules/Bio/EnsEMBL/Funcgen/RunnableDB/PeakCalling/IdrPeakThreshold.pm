
=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::IdrPeakThreshold

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::IdrPeakThreshold;

use warnings;
use strict;
use base 'Bio::EnsEMBL::Hive::Process';

use Data::Dumper;

sub run {
  my $self = shift;
  my $species        = $self->param_required('species');
  my $execution_plan = $self->param_required('execution_plan');
  my $idr_result     = $self->param_required('idr_result');

  my $experiment_name = $execution_plan->{idr}->{name};
  my $idr_strategy    = $execution_plan->{idr}->{strategy};
  
  my $experiment_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'experiment');
  my $experiment = $experiment_adaptor->fetch_by_name($experiment_name);
  
  if (! defined $experiment) {
    die;
  }
  
  my @pairwise_idr_values = map { $_->{idr_num_peak_threshold} } @$idr_result;
  
  my @failed_idr_pars;
  IDR_PAIR:
  foreach my $current_idr_result (@$idr_result) {
    if (! $current_idr_result->{idr_failed}) {
      next IDR_PAIR
    };
    my $peak_calling_pair = $current_idr_result->{peak_calling_pair};
    my $error_message     = $current_idr_result->{error_message};
    
    my $first_alignment_name  = $peak_calling_pair->[0]->{alignment_name};
    my $second_alignment_name = $peak_calling_pair->[1]->{alignment_name};
    
    push @failed_idr_pars, "($first_alignment_name, $second_alignment_name, $error_message)";
  }
  
  use List::Util qw( max );
  my $idr_num_peak_threshold = max ( @pairwise_idr_values );
  
  use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw (
    RUN_IDR_ON_BIOLOGICAL_REPLICATES 
    RUN_IDR_ON_TECHNICAL_REPLICATES  
  );
  use Bio::EnsEMBL::Funcgen::Idr qw ( 
    IDR_ON_BIOLOGICAL_REPLICATES
    IDR_ON_TECHNICAL_REPLICATES
  );

  my $idr_type;
  if ($idr_strategy eq RUN_IDR_ON_BIOLOGICAL_REPLICATES) {
    $idr_type = IDR_ON_BIOLOGICAL_REPLICATES;
  }
  if ($idr_strategy eq RUN_IDR_ON_TECHNICAL_REPLICATES) {
    $idr_type = IDR_ON_TECHNICAL_REPLICATES;
  }
  if (! defined $idr_type) {
    die;
  }
  
  my $failed_idr_pair_string = undef;
  
  if (@failed_idr_pars) {
    $failed_idr_pair_string = join "\n", @failed_idr_pars;
  }
  
  my $idr = Bio::EnsEMBL::Funcgen::Idr->new(
    -experiment_id    => $experiment->dbID,
    -max_peaks        => $idr_num_peak_threshold,
    -type             => $idr_type,
    -failed_idr_pairs => $failed_idr_pair_string,
  );
  my $idr_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'Idr');
  $idr_adaptor->store($idr);
  return;
}

1;
