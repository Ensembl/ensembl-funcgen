=head1 NAME

Bio::EnsEMBL::Funcgen::Hive::AddJobParameter

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::Hive::AddJobParameter;

use base ('Bio::EnsEMBL::Hive::Process');
use Data::Dumper;
use warnings;
use strict;

my $pipeline_branch = 2;

sub run {
  my $self = shift;
  
  my $input_job = $self->input_job;
  $input_job->autoflow(0);

  my $add = $self->param('add');
#   my %add_parameters = split ',', $add;
  
  use Bio::EnsEMBL::Hive::Utils qw( destringify );
  my $input_id = destringify($input_job->input_id);
  my %new_input_id = (%$input_id, %$add);
  
  $self->dataflow_output_id(\%new_input_id, $pipeline_branch);
  return;
}

1;
