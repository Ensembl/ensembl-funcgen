package Bio::EnsEMBL::Funcgen::Hive::CreateJobBatchUsingNewGroupingMechanism;

use base Bio::EnsEMBL::Hive::Process;
use strict;

sub _fetch_names_of_unprocessed_experiments {

  my $self = shift;
  my $funcgen_dbc = shift;
  
  my $helper =
    Bio::EnsEMBL::Utils::SqlHelper->new( -DB_CONNECTION => $funcgen_dbc );
    
  my $sql = '
  select 
    experiment.name as experiment_name
  from 
    experiment 
    join feature_type using (feature_type_id) 
    join input_subset using (experiment_id) 
    join input_subset_tracking using (input_subset_id) 
  where 
    experiment.is_control!=1
  and local_url not like "/lustre%"
  and experiment.experiment_id not in (select distinct experiment_id from feature_set where experiment_id is not null)
  order by
    control_id
  ;';

  my $arr_ref = $helper->execute(
    -SQL => $sql,
    -CALLBACK => sub {
      my $row = shift;
      return $row->[0];
    },
  );
  return $arr_ref;
}

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
  
  # An array reference to experiment names.
  my $experiment_name = $self->param('experiment_name');
  my $species         = $self->param('species');
  
  if (ref $experiment_name ne 'ARRAY') {

    $experiment_name = $self->_fetch_names_of_unprocessed_experiments(
      Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen')->dbc
    );

  }
#   print Dumper(ref $experiment_name ne 'ARRAY');
#   print Dumper($experiment_name);
#   die;

  use Bio::EnsEMBL::Utils::Logger;
  my $logger = Bio::EnsEMBL::Utils::Logger->new();
  $logger->init_log;

  use Bio::EnsEMBL::Registry;
  use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
  
  my $experiment_adaptor   = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'Experiment');
  my $input_subset_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'InputSubSet');

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

  my %control_name_to_signal_name = build_control_name_to_signal_name_hash($experiment_adaptor, $experiment_name);
  
  use Data::Dumper;
  
#   $logger->warning("\n-----------------------------------------------------------------\n");
#   $logger->warning(Dumper(\%control_name_to_signal_name));
#   $logger->warning("\n-----------------------------------------------------------------\n");
  
  if (values %control_name_to_signal_name == 0) {
    $logger->warning("No signal experiments found for\n" . Dumper($experiment_name));
  }
  
  my $number_of_experiment_groups_seeded = 0;
#   my $max_number_of_experiment_groups = 5;
  
  # Iterator over all control experiment names, for every signal experiment, 
  # create an ArrayRef of input_subset_ids that belong to this experiment.
  #
  # Add the control experiment's input_subet_id as an extra hash key.
  #
  # This is the format that DefineResultSets expects.
  # 
  BATCH: foreach my $current_control (keys %control_name_to_signal_name) {

    my $signal = $control_name_to_signal_name{$current_control};
    my $batch_job_definition = {};
    
    foreach my $current_signal_name (@$signal) {    
      $batch_job_definition->{$current_signal_name} = $fetch_input_subset_ids_by_experiment_name->($current_signal_name);
    }
    $batch_job_definition->{'controls'} = $fetch_input_subset_ids_by_experiment_name->($current_control);
    
#     print Dumper($batch_job_definition);
#     die;
    
    $self->dataflow_output_id({
      input_subset_ids => $batch_job_definition
    }, 2);
    
    $number_of_experiment_groups_seeded++;
    
#     if ($number_of_experiment_groups_seeded == $max_number_of_experiment_groups) {
#       last BATCH;
#     }
  }
  
  $logger->finish_log;
  
  return;
}

=head2 build_control_name_to_signal_name_hash

  Input: ArrayRef of signal experiment names
  Output: A hash mapping from control experiment names to ArrayRefs of 
    the signal experiment names.

=cut
sub build_control_name_to_signal_name_hash {

  my $experiment_adaptor = shift;
  my $experiment_name    = shift;
#   my $experiment_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'Experiment');

  my %control_name_to_signal_name;
  foreach my $current_experiment_name (@$experiment_name) {
    my $experiment = $experiment_adaptor->fetch_by_name($current_experiment_name);
    if (!defined $experiment) {
      use Carp;
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
