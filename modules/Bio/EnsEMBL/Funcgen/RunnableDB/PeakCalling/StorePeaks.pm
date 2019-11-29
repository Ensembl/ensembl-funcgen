package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::StorePeaks;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;
use Carp;
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw ( :all );
use JSON::XS qw(decode_json);
sub run {

  my $self = shift;
  my $plan                   = $self->param_required('execution_plan');
  my $species                = $self->param_required('species');
  my $peak_calling_succeeded = $self->param_required('peak_calling_succeeded');
  my $error_message          = $self->param('error_message');
  my $peaks_to_load_file     = $self->param('peaks_to_load_file');
  
  if ($peak_calling_succeeded && ! defined $peaks_to_load_file) {
    confess("peaks_to_load_file parameter is mandatory, if the peak calling succeeded!");
  }
  
  use Bio::EnsEMBL::Funcgen::PeakCallingPlan::ExecutionPlanUtils qw (
        lock_execution_plan
        resolve_nonterminal_symbols
  );
  my $plan_expanded = resolve_nonterminal_symbols($plan);
  lock_execution_plan($plan_expanded);
  print Dumper($plan_expanded);
  
  my $experiment_name   = $plan_expanded->{idr}->{name};
  my $name              = $plan_expanded->{call_peaks}->{name};
  my $display_label     = $plan_expanded->{call_peaks}->{display_label};
  my $feature_type_name = $plan_expanded->{call_peaks}->{feature_type};
  my $logic_name        = $plan_expanded->{call_peaks}->{analysis};
  my $signal_alignment_name = $plan_expanded
    ->{call_peaks}
    ->{input}
    ->{signal}
    ->{name}
  ;
  my $control_alignment_name;

  if ($plan_expanded->{meta_data}->{experiment_has_control} eq FALSE) {
    $control_alignment_name = undef;
  } else {
    $control_alignment_name = $plan_expanded
    ->{call_peaks}
    ->{input}
    ->{control}
    ->{name}
  ;
  }

  my $slice_adaptor        = Bio::EnsEMBL::Registry->get_adaptor($species, 'core',    'Slice');
  my $peak_adaptor         = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'Peak');
  my $peak_calling_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'PeakCalling');
  my $feature_type_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'FeatureType');
  my $analysis_adaptor     = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'Analysis');
  my $alignment_adaptor    = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'Alignment');
  my $experiment_adaptor   = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'Experiment');
  
  my $analysis          = $analysis_adaptor     ->fetch_by_logic_name($logic_name);
  my $feature_type      = $feature_type_adaptor ->fetch_by_name($feature_type_name);
  my $signal_alignment  = $alignment_adaptor    ->fetch_by_name($signal_alignment_name);
  my $experiment        = $experiment_adaptor   ->fetch_by_name($experiment_name);
  
  if (! defined $analysis) {
    die("Can't find analysis with name $logic_name!");
  }
  if (! defined $feature_type) {
    die("Can't find feature type with name $feature_type_name!");
  }
  if (! defined $signal_alignment) {
    die("Can't find alignment with name $signal_alignment_name!");
  }
  if (! defined $experiment) {
    die("Can't find experiment with name $experiment_name!");
  }
  
  my $control_alignment_id;
  if ($control_alignment_name) {
    my $control_alignment = $alignment_adaptor->fetch_by_name($control_alignment_name);
    
    if (! defined $control_alignment) {
        confess("Couldn't fetch control alignment with name $control_alignment_name for experiment $experiment_name!");
    }
    
    $control_alignment_id = $control_alignment->dbID,
  } else {
    $control_alignment_id = undef;
  }

  my $peak_calling = $peak_calling_adaptor->fetch_by_name($name);
  
  if (! defined $peak_calling) {
  
    use Bio::EnsEMBL::Funcgen::PeakCalling;
    $peak_calling = Bio::EnsEMBL::Funcgen::PeakCalling->new(
      -name                 => $name,
      -display_label        => $display_label,
      -feature_type_id      => $feature_type->dbID,
      -analysis_id          => $analysis->dbID,
      -signal_alignment_id  => $signal_alignment->dbID,
      -control_alignment_id => $control_alignment_id,
      -experiment_id        => $experiment->dbID,
      -epigenome_id         => $experiment->epigenome->dbID,
    );
    $peak_calling_adaptor->store($peak_calling);
  }
  
  if (!$peak_calling_succeeded) {
    
    $peak_calling->run_failed(! $peak_calling_succeeded);
    $peak_calling->error_message($error_message);
    $peak_calling_adaptor->update($peak_calling);
    return;
    
  }

  open my $peak_in, '<' . $peaks_to_load_file || die("Couldn't open file $peaks_to_load_file for reading!");
  while (my $line = <$peak_in>){
      chomp $line;
      my %peak_to_store = %{ decode_json $line };
      store_peaks_in_db($self, $slice_adaptor, $peak_adaptor, \%peak_to_store, $peak_calling);
  }
  $peak_in->close;

  return;
}

sub store_peaks_in_db {
  my $self = shift;
  my $slice_adaptor = shift;
  my $peak_adaptor = shift;
  my $peak_hash = shift;
  my $peak_calling = shift;
    
  use Hash::Util qw( lock_hash );
  lock_hash(%$peak_hash);
  
  if (! exists  $peak_hash->{-seq_region}) {
    warn("seq_region wasn't defined for this peak. Skipping.");
    return;
  }
  
  my $slice = $slice_adaptor->fetch_by_region(
    'toplevel', 
    $peak_hash->{-seq_region}
  );
  if (! defined $slice) {
    warn("Couldn't fetch slice: " . $peak_hash->{-seq_region});
    return
  }
  use Bio::EnsEMBL::Funcgen::Peak;
  my $peak = Bio::EnsEMBL::Funcgen::Peak->new(
    -slice             => $slice,
    -seq_region_start  => $peak_hash->{-start},
    -seq_region_end    => $peak_hash->{-end},
    -seq_region_strand => $peak_hash->{-strand},
    -summit            => $peak_hash->{-summit},
    -score             => $peak_hash->{-score},
    -peak_calling      => $peak_calling,
  );
  
  eval {
    $peak_adaptor->store($peak);
  };
  if ($@) {
    $self->throw(
      $@ . "\n\n"
      . Dumper($peak_hash)
      . Dumper($peak)
    );
  }
  # See, if this fixes the rare occasion of memory leaks.
  undef($peak);
  return;
}

1;
