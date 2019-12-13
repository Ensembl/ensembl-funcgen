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

=head1 SYNOPSIS

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::Report::binnedValueReport;

use strict;
use warnings;

use Role::Tiny;

sub _generate_bins {

  my $self = shift;
  
  my $param = shift;
  
  my $min      = $param->{min};
  my $max      = $param->{max};
  my $num_bins = $param->{num_bins};

  my @bins
    =
      map { $_ / $num_bins * ($max - $min) + $min }
      0..$num_bins
    ;
  return \@bins
}

sub _grep {

  my $self = shift;
  
  my $object_list     = shift;
  my $filter_callback = shift;
  
  if (ref $filter_callback ne 'CODE') {
    use Carp;
    confess("filter_callback must be a code reference!");
  }
  
  my $filtered_object_list = [ 
    grep { $filter_callback->($_) } @$object_list 
  ];
  
  return $filtered_object_list;
}

sub _merge_hashrefs {

  my $self = shift;
  
  my $hashref1 = shift;
  my $hashref2 = shift;
  
  my $merged = {
    %$hashref1,
    %$hashref2,
  };
  return $merged;
}

sub _compute_dataset_with_bin_counts {
  
  my $self  = shift;
  my $param = shift;
  
  my $object_list             = $param->{object_list};
  my $static_values           = $param->{static_values};
  my $filter_callback         = $param->{object_filter};
  my $compute_values_sub_name = $param->{compute_value_method};
  
  if (ref $object_list ne 'ARRAY') {
    confess("object_list parameter must be an array reference!");
  }

  if (@$object_list == 0) {
    warn("Objects list is empty!");
  }

  if (! defined $filter_callback) {
    $filter_callback = sub { return 1; };
  }

  if (! defined $compute_values_sub_name) {
    $compute_values_sub_name = '_compute_values_from_object_list';
  }

  my $filtered_object_list = $self->_grep($object_list, $filter_callback);
  
  if (@$filtered_object_list == 0) {
    warn("No objects left after applying filter!");
  }
  
  my $dynamic_values
    = $self->$compute_values_sub_name($filtered_object_list);

  my $final_dataset = $self->_merge_hashrefs(
    $static_values,
    $dynamic_values,
  );
  return $final_dataset;
}

sub _count_values_per_bin {

  my $self   = shift;
  my $bins   = shift;
  my $values = shift;

  use Statistics::Descriptive;
  my $stat = Statistics::Descriptive::Full->new();

  $stat->add_data( @$values );
  my $f = $stat->frequency_distribution_ref( $bins );

  my @counts;
  for (sort {$a <=> $b} keys %$f) {
    push @counts, $f->{$_};
  }
  return \@counts;
}

1;
