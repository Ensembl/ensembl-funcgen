=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

package Bio::EnsEMBL::Funcgen::DBSQL::GenericFeatureAdaptor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw );
use base 'Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor';

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::DBSQL::GenericAdaptorMethods';

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  $self->init_generic_adaptor(@args);
  $self->init(@args);
  return $self;
}

sub object_class {
  my $self = shift;
  throw('object_class has to be overwritten by the subclass in ('. ref $self .')!');
}

use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_get_or_set
);

sub slice_cache {
  my $self  = shift;
  my $value = shift;
  return $self->_generic_get_or_set('slice_cache', $value);
}

sub init {
  my $self = shift;
  $self->slice_cache({});
}

sub _objs_from_sth {
    my ($self, $sth, $mapper, $dest_slice) = @_;
    
    my @features;
    FEATURE: while ( my $row = $sth->fetchrow_hashref ) {
        my $feature = $self->objectify($row);
        $self->_load_dependencies($feature, $dest_slice, $mapper);
        push @features, $feature;
    }
    return \@features;
}

sub _load_dependencies {
    my $self       = shift;
    
    my $feature    = shift;
    my $dest_slice = shift;
    my $mapper     = shift;
    
    (
      my $feature_slice, 
      my $feature_start, 
      my $feature_end, 
      my $feature_strand
    ) = $self->_create_feature_slice($feature, $dest_slice, $mapper);
    
    $feature->slice  ($feature_slice);
    $feature->strand ($feature_strand);
    $feature->start  ($feature_start);
    $feature->end    ($feature_end);

    return $feature;
}

sub _create_feature_slice {
  my ($self, $feature, $dest_slice, $mapper) = @_;
  
  my $seq_region_id    = $feature->seq_region_id;
  my $strand           = $feature->seq_region_strand;
  my $seq_region_start = $feature->seq_region_start;
  my $seq_region_end   = $feature->seq_region_end;
  
  if ($mapper) {
    throw("This is a data issue. There should be no features on non-toplevel seq regions.");
  }

  # If a destination slice was provided convert the coords
  # If the destination slice starts at 1 and is forward strand, nothing needs doing
  if ($dest_slice) {
  
    my $dest_slice_start   = $dest_slice->start;
    my $dest_slice_end     = $dest_slice->end;
    my $dest_slice_strand  = $dest_slice->strand;

    unless ($dest_slice_start == 1 && $dest_slice_strand == 1) {
      if ($dest_slice_strand == 1) {
        $seq_region_start = $seq_region_start - $dest_slice_start + 1;
        $seq_region_end   = $seq_region_end   - $dest_slice_start + 1;
      } else {
        my $tmp_seq_region_start = $seq_region_start;
        $seq_region_start        = $dest_slice_end - $seq_region_end       + 1;
        $seq_region_end          = $dest_slice_end - $tmp_seq_region_start + 1;
        $strand      *= -1;
      }
    }
    return $dest_slice, $seq_region_start, $seq_region_end, $strand;
  }

  my $slice_adaptor    = $self->db->dnadb->get_SliceAdaptor;
  my $slice = $self->slice_cache->{$seq_region_id};

  if (! $slice) {
    $slice = $slice_adaptor->fetch_by_seq_region_id($seq_region_id);
    $self->slice_cache->{$seq_region_id} = $slice;
  }
  return $slice, $seq_region_start, $seq_region_end, $strand;
}

1;
