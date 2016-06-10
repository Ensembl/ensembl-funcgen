package Bio::EnsEMBL::Funcgen::Hive::QcFastQcInputIdsFromInputSet;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';

sub run {
  my $self = shift;
  my $input_subset_ids = $self->param('input_subset_ids');
  
  die ('Type error') unless(ref $input_subset_ids eq 'HASH');
  
  my @keys = keys %$input_subset_ids;
  
  my @input_subset_id;
  foreach my $current_key (@keys) {
    my $current_input_subset_id = $input_subset_ids->{$current_key};

    die ('Type error') unless(ref $current_input_subset_id eq 'ARRAY');
    push @input_subset_id, @$current_input_subset_id;
  }
  
  foreach my $current_input_subset_id (@input_subset_id) {
    $self->dataflow_output_id(
      { input_subset_id => $current_input_subset_id }, 
      2
    );
  }
  return;
}

1;
