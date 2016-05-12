=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Funcgen::RegulatoryActivity

=head1 SYNOPSIS
=head1 DESCRIPTION
=head1 SEE ALSO
=cut

package Bio::EnsEMBL::Funcgen::RegulatoryActivity;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw deprecate );

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self =  bless({}, $class);
  return $self;
}

=head2 feature_set_id
=cut
sub feature_set_id {
  my $self = shift;
  my $feature_set_id = shift;

  if(defined $feature_set_id) {
    $self->{'_feature_set_id'}  = $feature_set_id;
  }
  return $self->{'_feature_set_id'};
}

=head2 activity
=cut
sub activity {
  my $self = shift;
  my $activity = shift;

  if(defined $activity) {
    $self->{'_activity'}  = $activity;
  }
  return $self->{'_activity'};
}

=head2 set_regulatory_evidence
=cut
sub set_regulatory_evidence {
  my $self = shift;
  my $regulatory_evidence = shift;

  $self->{'_regulatory_evidence'}  = $regulatory_evidence;
  return;
}

sub db {
  my $self = shift;
  my $db   = shift;

  if ($db) {
    $self->{'_db'}  = $db;
  }
  return $self->{'_db'};
}

=head2 get_regulatory_evidence
=cut
sub get_regulatory_evidence {
  my $self = shift;
  my $feature_class = shift;

  if (! defined $feature_class) {
    throw("Missing parameter for feature_class!");
  }
  if ($feature_class eq 'annotated') {
    return $self->supporting_annotated_features;
  }
  if ($feature_class eq 'motif') {
    return $self->supporting_motif_features;
  }
  throw("Unknown feature_class $feature_class!");
}

sub supporting_annotated_features {
  my $self = shift;
  
  my @id = $self->supporting_annotated_feature_ids;
  my @annotated_feature = map {
    $self->db->get_AnnotatedFeatureAdaptor->fetch_by_dbID($_);
  } @id;
  
  return \@annotated_feature;
}

sub supporting_motif_features {
  my $self = shift;
  
  my @id = $self->supporting_motif_feature_ids;
  my @motif_feature = map {
    $self->db->get_MotifFeatureAdaptor->fetch_by_dbID($_);
  } @id;
  
  return \@motif_feature;
}

sub supporting_annotated_feature_ids {
  my $self = shift;
  return keys %{$self->{'_regulatory_evidence'}->{'annotated'}};
}

sub supporting_motif_feature_ids {
  my $self = shift;
  return keys %{$self->{'_regulatory_evidence'}->{'motif'}};
}

sub get_underlying_structure {
  my $self = shift;
  
  my @motif_feature_loci;
  foreach my $current_motif_feature (@{$self->supporting_motif_features}) {
    push @motif_feature_loci, (
      0 + $current_motif_feature->start, 
      0 + $current_motif_feature->end
    );
  }

  return \@motif_feature_loci;
}


1;


