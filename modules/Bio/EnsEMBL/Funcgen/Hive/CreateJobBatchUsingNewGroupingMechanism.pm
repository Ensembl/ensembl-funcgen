package Bio::EnsEMBL::Funcgen::Hive::CreateJobBatchUsingNewGroupingMechanism;

use base Bio::EnsEMBL::Hive::Process;
use strict;

=head2

Builds something like this:

{
  'K562:hist:BR1_H3K27ac_3526' => [3306,3332,3333,3339,3380,3402],
  'K562:hist:BR1_H3K27me3_3526' => [3219,3248,3447],
  'K562:hist:BR1_H3K4me3_3526' => [3233,3301,3366,3403,3408,3444],
  'controls' => [3461]
}

=cut

sub run {
  my $self = shift; 
  
  my $out_db_url      = $self->param('out_db_url');
  my $experiment_name = $self->param('experiment_name');
  my $species         = $self->param('species');

  use Bio::EnsEMBL::Registry;
  use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
  Bio::EnsEMBL::Registry->load_registry_from_url(
    "${out_db_url}?group=funcgen&species=$species", 
    1
  );
  
  my $experiment_adaptor   = Bio::EnsEMBL::Registry->get_adaptor('homo_sapiens', 'funcgen', 'Experiment');
  my $input_subset_adaptor = Bio::EnsEMBL::Registry->get_adaptor('homo_sapiens', 'funcgen', 'InputSubSet');

  # For a given experiment name, creates a list of all input_subset_ids 
  # belonging to this experiment.
  #
  my $fetch_input_subset_ids_by_experiment_name = sub {

      my $experiment_name = shift;
      
      my $input_subset = $input_subset_adaptor->fetch_all_by_Experiments(
	$experiment_adaptor->fetch_by_name($experiment_name)
      );
      my @input_subset_id = map { $_->dbID } @$input_subset;
      return \@input_subset_id;
  };

  my %control_name_to_signal_name = build_control_name_to_signal_name_hash($experiment_name);
  
  # Iterator over all control experiment names, for every signal experiment, 
  # create an ArrayRef of input_subset_ids that belong to this experiment.
  #
  # Add the control experiment's input_subet_id as an extra hash key.
  #
  # This is the format that DefineResultSets expects.
  # 
  foreach my $current_control (keys %control_name_to_signal_name) {

    my $signal = $control_name_to_signal_name{$current_control};
    my $batch_job_definition = {};
    
    foreach my $current_signal_name (@$signal) {    
      $batch_job_definition->{$current_signal_name} = $fetch_input_subset_ids_by_experiment_name->($current_signal_name);
    }
    $batch_job_definition->{'controls'} = $fetch_input_subset_ids_by_experiment_name->($current_control);
    
    $self->dataflow_output_id({
      input_subset_ids => $batch_job_definition
    }, 2);
  }
  return;
}

=head2 build_control_name_to_signal_name_hash

  Input: ArrayRef of signal experiment names
  Output: A hash mapping from control experiment names to ArrayRefs of 
    the signal experiment names.

=cut
sub build_control_name_to_signal_name_hash {

  my $experiment_name = shift;
  my $experiment_adaptor = Bio::EnsEMBL::Registry->get_adaptor('homo_sapiens', 'funcgen', 'Experiment');

  my %control_name_to_signal_name;
  foreach my $current_experiment_name (@$experiment_name) {
    my $experiment = $experiment_adaptor->fetch_by_name($current_experiment_name);
    if (!defined $experiment) {
      confess("Can't find experiment with name $current_experiment_name!");
    }
    if (! $experiment->is_control) {
      my $control_experiment = $experiment->get_control;
      if (! exists $control_name_to_signal_name{$control_experiment->name}) {
	$control_name_to_signal_name{$control_experiment->name} = [];
      }
      push @{$control_name_to_signal_name{$control_experiment->name}}, $experiment->name;
    }
  }
  return %control_name_to_signal_name;
}

1;
