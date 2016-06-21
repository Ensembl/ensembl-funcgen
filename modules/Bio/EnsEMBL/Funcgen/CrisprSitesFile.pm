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

Bio::EnsEMBL::Funcgen::CrisprSitesFile

=head1 SYNOPSIS
=head1 DESCRIPTION
=cut

package Bio::EnsEMBL::Funcgen::CrisprSitesFile;

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
    file
    file_type
    analysis
    adaptor
  );
  
  my (
    $dbID,
    $name,
    $file,
    $file_type,
    $analysis,
    $adaptor,
  )
    = rearrange([ @field ], @_);

  $self->dbID            ($dbID);
  $self->name            ($name);
  $self->file            ($file);
  $self->file_type       ($file_type);
  $self->_analysis       ($analysis);
  $self->adaptor         ($adaptor);

  return $self;
}

sub dbID           { return shift->_generic_get_or_set('dbID',            @_) }
sub name           { return shift->_generic_get_or_set('name',            @_) }
sub file           { return shift->_generic_get_or_set('file',            @_) }
sub file_type      { return shift->_generic_get_or_set('file_type',       @_) }
sub _analysis      { return shift->_generic_get_or_set('_analysis',       @_) }
sub adaptor        { return shift->_generic_get_or_set('adaptor',         @_) }

sub get_Analysis {
  my $self = shift;
  return $self->_analysis
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


