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

sub db {
  my $self = shift;
  my $db   = shift;

  if ($db) {
    $self->{'_db'}  = $db;
  }
  return $self->{'_db'};
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

=head2 regulatory_evidence
=cut
sub regulatory_evidence {
  my $self = shift;
  my $regulatory_evidence = shift;

  if ($regulatory_evidence) {
    $self->{'_regulatory_evidence'} = $regulatory_evidence;
  }
  if (! $self->{'_regulatory_evidence'}) {
    use Bio::EnsEMBL::Funcgen::RegulatoryEvidence;
    my $regulatory_evidence = Bio::EnsEMBL::Funcgen::RegulatoryEvidence->new;
    $regulatory_evidence->db($self->db);
    $self->{'_regulatory_evidence'} = $regulatory_evidence;
  }
  return $self->{'_regulatory_evidence'};
}

1;
