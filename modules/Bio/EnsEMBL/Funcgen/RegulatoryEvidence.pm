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

Bio::EnsEMBL::Funcgen::RegulatoryEvidence

=head1 SYNOPSIS
=head1 DESCRIPTION
=head1 SEE ALSO
=cut

package Bio::EnsEMBL::Funcgen::RegulatoryEvidence;

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

=head2 set_regulatory_evidence
=cut
sub set_regulatory_evidence {
  my $self = shift;
  my $regulatory_evidence = shift;

#   $self->{'_regulatory_evidence'}  = $regulatory_evidence;
  
  $self->{'_annotated'}  = $regulatory_evidence->{'annotated'};
  $self->{'_motif'}      = $regulatory_evidence->{'motif'};
  
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

sub supporting_annotated_features {
  my $self = shift;
  
  # Stopgap fix, web is getting this error:
  #
  # MSG: Trying to assign a non numerical parameter to an integer value in the database
  # STACK Bio::EnsEMBL::DBSQL::BaseAdaptor::bind_param_generic_fetch /ensemblweb/ap5/www/reg_fixes/ensembl/modules/Bio/EnsEMBL/DBSQL/BaseAdaptor.pm:381
  # STACK Bio::EnsEMBL::DBSQL::BaseAdaptor::_uncached_fetch_by_dbID /ensemblweb/ap5/www/reg_fixes/ensembl/modules/Bio/EnsEMBL/DBSQL/BaseAdaptor.pm:670
  # STACK Bio::EnsEMBL::DBSQL::BaseAdaptor::fetch_by_dbID /ensemblweb/ap5/www/reg_fixes/ensembl/modules/Bio/EnsEMBL/DBSQL/BaseAdaptor.pm:653
  # STACK Bio::EnsEMBL::Funcgen::RegulatoryEvidence::supporting_annotated_features /ensemblweb/ap5/www/reg_fixes/ensembl-funcgen/modules/Bio/EnsEMBL/Funcgen/RegulatoryEvidence.pm:76
  #
  return;
  
  my @id = $self->supporting_annotated_feature_ids;
  my @annotated_feature = map {
    $self->db->get_AnnotatedFeatureAdaptor->fetch_by_dbID($_);
  } @id;
  
  return \@annotated_feature;
}

sub supporting_motif_features {
  my $self = shift;
  
  # Stopgap fix, web is getting the same error here as above with 
  # supporting_annotated_features:
  #
  return;
  my @id = $self->supporting_motif_feature_ids;
  my @motif_feature = map {
    $self->db->get_MotifFeatureAdaptor->fetch_by_dbID($_);
  } @id;
  
  return \@motif_feature;
}

sub add_supporting_annotated_feature_id {
  my $self = shift;
  my $annotated_feature_id = shift;
  
  if (! defined $annotated_feature_id) {
    return;
  }
  if (ref $annotated_feature_id eq 'ARRAY') {
    foreach my $id (@$annotated_feature_id) {
      $self->add_supporting_annotated_feature_id($id);
    }
  }
  
  $self->{'_annotated'}->{0 + $annotated_feature_id} = undef;
}

sub add_supporting_motif_feature_id {
  my $self = shift;
  my $motif_feature_id = shift;
  
  if (! defined $motif_feature_id) {
    return;
  }
  if (ref $motif_feature_id eq 'ARRAY') {
    foreach my $id (@$motif_feature_id) {
      $self->add_supporting_motif_feature_id($id);
    }
  }
  
  $self->{'_motif'}->{0 + $motif_feature_id} = undef;
}

sub supporting_annotated_feature_ids {
  my $self = shift;
  if (! defined $self->{'_annotated'}) {
	return 
  }
  my @supporting_annotated_feature_ids = keys %{$self->{'_annotated'}};
  return \@supporting_annotated_feature_ids;
}

sub supporting_motif_feature_ids {
  my $self = shift;
  if (! defined $self->{'_motif'}) {
	return;
  }
  my @supporting_motif_feature_ids = keys %{$self->{'_motif'}};
  return \@supporting_motif_feature_ids;
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
