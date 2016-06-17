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

Bio::EnsEMBL::Funcgen::RegulatoryBuild

=head1 SYNOPSIS
=head1 DESCRIPTION
=cut

package Bio::EnsEMBL::Funcgen::RegulatoryBuild;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw deprecate );

=head2 new
=cut
sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = bless {}, $class;

  my @field = qw(
    dbID
    name
    version
    initial_release_date
    last_annotation_update
    feature_type_id
    analysis_id
    is_current
    epigenomes
  );
  
  my (
    $dbID,
    $name,
    $version,
    $initial_release_date,
    $last_annotation_update,
    $feature_type_id,
    $analysis_id,
    $is_current,
    $epigenomes
  )
    = rearrange([ @field ], @_);

  $self->dbID($dbID);
  $self->name($name);
  $self->version($version);
  $self->initial_release_date($initial_release_date);
  $self->last_annotation_update($last_annotation_update);
  $self->feature_type_id($feature_type_id);
  $self->analysis_id($analysis_id);
  $self->is_current($is_current);
  $self->epigenomes($epigenomes);

  return $self;
}

sub dbID                   { return shift->_generic_get_or_set('dbID',                   @_) }
sub name                   { return shift->_generic_get_or_set('name',                   @_) }
sub version                { return shift->_generic_get_or_set('version',                @_) }
sub initial_release_date   { return shift->_generic_get_or_set('initial_release_date',   @_) }
sub last_annotation_update { return shift->_generic_get_or_set('last_annotation_update', @_) }
sub feature_type_id        { return shift->_generic_get_or_set('feature_type_id',        @_) }
sub analysis_id            { return shift->_generic_get_or_set('analysis_id',            @_) }
sub is_current             { return shift->_generic_get_or_set('is_current',             @_) }
sub epigenomes             { return shift->_generic_get_or_set('epigenomes',             @_) }

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


