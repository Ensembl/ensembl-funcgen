package Bio::EnsEMBL::Funcgen::PeakCallingPlan::ExecutionPlanUtils;

use warnings;
use strict;
use Carp;

use base qw( Exporter );
use vars qw( @EXPORT_OK );
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw ( :all );

@EXPORT_OK = qw(
  create_ref
  resolve_nonterminal_symbols
  lock_execution_plan
  summarise
);

sub summarise {

  my $execution_plan = shift;
  
  use Data::Dumper;
  my $d = Data::Dumper
    ->new(
      [ $execution_plan ], 
      ['e']
    )
    ->Deepcopy(1)
    ->Sortkeys(1);
  return $d->Dump;
}

sub create_ref {
  my $plan = shift;
  
  my $name = $plan->{name};
  my $type = $plan->{type};
  
  if (! defined $name) {
    die;
  }
  if (! defined $type) {
    die;
  }
  
  my $ref = '%' . $type . ':' . $name . '%';
  
  # For each experiment, there is only one of these, so there is no
  # need to name them to keep them separate.
  #
  if ($type eq IDR_ANALYSIS || $type eq CALL_PEAKS_ANALYSIS) {
    $ref = '%' . $type . ':' . 'NA' . '%';
  }
  return $ref;
}

sub lock_execution_plan {

  my $execution_plan = shift;
  
  if (! ref $execution_plan) {
    return;
  }
  
  if (ref $execution_plan eq 'HASH') {
  
    use Hash::Util qw ( lock_hash );
    lock_hash( %$execution_plan );
  
    my @keys = keys %$execution_plan;
    foreach my $current_key (@keys) {
      lock_execution_plan(
        $execution_plan->{$current_key}
      );
    }
  }
  
  if (ref $execution_plan eq 'ARRAY') {
  
    for(my $i=0; $i<@$execution_plan; $i++) {
      lock_execution_plan(
        $execution_plan->[$i]
      );
    }
  }
  return;
}

sub resolve_nonterminal_symbols {

  my $execution_plan = shift;
  
  use Storable qw(dclone);
  my $execution_plan_clone;
  
  eval {
    $execution_plan_clone = dclone($execution_plan);
  };
  if ($@) {
    confess(
      Dumper($execution_plan) . "\n"
      . $@
    );
  }

  _resolve_nonterminal_symbols(
    $execution_plan_clone, 
    $execution_plan_clone
  );
  return $execution_plan_clone;
}

sub _resolve_nonterminal_symbols {

  my $execution_plan     = shift;
  my $current_sub_branch = shift;
  
  if (ref $current_sub_branch eq 'HASH') {
  
    my @keys = keys %$current_sub_branch;
    foreach my $key (@keys) {
      $current_sub_branch->{$key} 
        = _resolve_nonterminal_symbols(
          $execution_plan, 
          $current_sub_branch->{$key},
        );
    }
  }

  if (ref $current_sub_branch eq 'ARRAY') {
  
    for(my $i=0; $i<@$current_sub_branch; $i++) {
      $current_sub_branch->[$i]
        = _resolve_nonterminal_symbols(
          $execution_plan, 
          $current_sub_branch->[$i]
        );
      }
  }

  return _resolve_nonterminal_symbols_scalar(
    $execution_plan, 
    $current_sub_branch
  );
}

sub _resolve_nonterminal_symbols_scalar {

  my $execution_plan = shift;
  my $current_value  = shift;
  
  if (! defined $current_value) {
    return;
  }
  
  (
    my $is_nonterminal,
    my $type,
    my $name,
  ) = check_is_nonterminal_symbol($current_value);
  
  if ($is_nonterminal) {
  
    my $product;
    if ($type eq IDR_ANALYSIS || $type eq CALL_PEAKS_ANALYSIS) {
        $product = $execution_plan->{$type};
    } else {
        $product = $execution_plan->{$type}->{$name};
    }
    
    if (! defined $product) {
      die("Can't resolve: $current_value - name = $name, type = $type");
    }
    _resolve_nonterminal_symbols(
      $execution_plan, 
      $product
    );
    return $product;
  }
  return $current_value;
}

sub check_is_nonterminal_symbol {
  my $symbol = shift;
  
  # Character class \w doesn't work, because there are experiments with a -
  # in them.
  #
  #my $is_nonterminal = $symbol =~ /^#(\w+):(\w+)#$/;
  my $is_nonterminal = $symbol =~ /^%([a-zA-Z0-9_-]+):([a-zA-Z0-9_-]+)%$/;
  
  my $type;
  my $name;
  
  if ($is_nonterminal) {
    $type = $1;
    $name = $2;
  }
  return $is_nonterminal, $type, $name;
}

1;
