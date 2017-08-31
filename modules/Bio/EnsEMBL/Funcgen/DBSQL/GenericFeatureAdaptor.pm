=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

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
use Bio::EnsEMBL::Utils::Exception qw( throw warning );

use base (
  'Bio::EnsEMBL::Funcgen::DBSQL::GenericAdaptor',
  'Bio::EnsEMBL::DBSQL::BaseFeatureAdaptor',
);

sub init {
  my $self = shift;
  $self->SUPER::init($@);
  $self->slice_cache({});
}

sub _objs_from_sth {
    my ($self, $sth, $mapper, $dest_slice) = @_;
    
    my @features;
    FEATURE: while ( my $row = $sth->fetchrow_hashref ) {

        my $object = $self->objectify($row);
        
        (
          my $feature_slice, 
          my $feature_start, 
          my $feature_end, 
          my $feature_strand
        ) = $self->_create_feature_slice($object, $dest_slice, $mapper);
        
        $object->slice  ($feature_slice);
        $object->strand ($feature_strand);
        $object->start  ($feature_start);
        $object->end    ($feature_end);
        
        push @features, $object;
    }
    return \@features;
}

sub _simple_accessors {
  my $self = shift;
  my $super_accessors = $self->SUPER::_simple_accessors;
  
  my $accessors = [
    @$super_accessors,
    { method_name => 'slice_cache', hash_key => '_slice_cache',  },
  ];
  return $accessors;
}

sub _create_feature_slice {
  my ($self, $feature, $dest_slice, $mapper) = @_;
  
  my $seq_region_id    = $feature->seq_region_id;
  my $strand           = $feature->strand;
  my $seq_region_start = $feature->seq_region_start;
  my $seq_region_end   = $feature->seq_region_end;
  
  my $slice_adaptor    = $self->db->dnadb->get_SliceAdaptor;

  my (
    $dest_slice_start,   $dest_slice_end, 
    $dest_slice_strand,
  );

  if ($dest_slice) {
    $dest_slice_start   = $dest_slice->start;
    $dest_slice_end     = $dest_slice->end;
    $dest_slice_strand  = $dest_slice->strand;
  }

  my $slice = $self->slice_cache->{$seq_region_id};

  if (! $slice) {
    $slice = $slice_adaptor->fetch_by_seq_region_id($seq_region_id);
    $self->slice_cache->{$seq_region_id} = $slice;
  }

  if ($mapper) {
    throw("This is a data issue. There should be no features on non-toplevel seq regions.");
  }

  # If a destination slice was provided convert the coords
  # If the destination slice starts at 1 and is forward strand, nothing needs doing
  if ($dest_slice) {

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
  }
  return $dest_slice, $seq_region_start, $seq_region_end, $strand;
}

1;
