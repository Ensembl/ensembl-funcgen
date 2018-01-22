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

  Bio::EnsEMBL::Funcgen::ProbeTranscriptMapping

=head1 SYNOPSIS
=head1 DESCRIPTION
=cut

package Bio::EnsEMBL::Funcgen::ProbeTranscriptMapping;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw deprecate );

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = bless {}, $class;

  my @field = qw(
    dbID
    probe_id
    stable_id
    description
    adaptor
  );
  
  my (
    $dbID,
    $probe_id,
    $stable_id,
    $description,
    $adaptor
  )
    = rearrange([ @field ], @_);

  $self->dbID            ($dbID);
  $self->stable_id       ($stable_id);
  $self->description     ($description);
  $self->probe_id        ($probe_id);
  $self->adaptor         ($adaptor);

  return $self;
}

sub dbID           { return shift->_generic_get_or_set('dbID',            @_) }
sub adaptor        { return shift->_generic_get_or_set('adaptor',         @_) }
sub stable_id      { return shift->_generic_get_or_set('stable_id',       @_) }
sub description    { return shift->_generic_get_or_set('description',     @_) }
sub probe_id       { return shift->_generic_get_or_set('probe_id',        @_) }

sub fetch_Probe {
  my $self = shift;
  my $probe = $self->adaptor->get_ProbeAdaptor->fetch_by_dbID($self->probe_id);
  return $probe;
}

sub _generic_get_or_set {
  my $self  = shift;
  my $name  = shift;
  my $value = shift;

  if(defined $value) {
    $self->{$name}  = $value;
  }
  return $self->{$name};
}

1;


