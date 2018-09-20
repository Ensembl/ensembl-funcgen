package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::SeedPairwise;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;

use constant {
  BRANCH_OUTPUT => 2,
};

sub run {

  my $self = shift;
  
  my $plan                    = $self->param_required('execution_plan');
  my $species                 = $self->param_required('species');
  my $permissive_peak_calling = $self->param_required('permissive_peak_calling');

  use Bio::EnsEMBL::Funcgen::PeakCallingPlan::ExecutionPlanUtils qw (
        lock_execution_plan
        resolve_nonterminal_symbols
  );
  my $plan_expanded = resolve_nonterminal_symbols($plan);
  lock_execution_plan($plan_expanded);
  
  my $pairs = $self->all_unique_pairs_no_self_pairing($permissive_peak_calling);
  
  $self->say_with_header("Got " . @$pairs . " pairs.", 1);
  
  foreach my $current_pair (@$pairs) {
  
    my $first_alignment  = $current_pair->[0]->{alignment_name};
    my $second_alignment = $current_pair->[1]->{alignment_name};
    
    my $first_file  = $current_pair->[0]->{peak_file};
    
    use File::Basename qw( dirname basename );
    my $tempdir = dirname( $first_file );

    my $idr_output_file = $tempdir . '/' . $first_alignment . '_vs_' . $second_alignment . '.idrPeaks.bed';

    use File::Path qw( make_path );
    make_path($tempdir);

    $self->dataflow_output_id( 
      {
        pair            => $current_pair,
        tempdir         => $tempdir,
        idr_output_file => $idr_output_file,
      }, 
      BRANCH_OUTPUT
    );
  }
}

sub all_unique_pairs_no_self_pairing {

  my $self  = shift;
  my $array = shift;

  my @pairs;
  
  for (my $i=0; $i<@$array; $i++) {
    for (my $j=$i+1; $j<@$array; $j++) {
      push @pairs, [
        $array->[$i],
        $array->[$j]
      ];
    }
  }
  return \@pairs;
}

1;
