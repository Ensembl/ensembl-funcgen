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

Bio::EnsEMBL::Funcgen::SegmentationFile

=head1 SYNOPSIS
=head1 DESCRIPTION
=cut

package Bio::EnsEMBL::Funcgen::SegmentationFile;

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
    logic_name
    description
    display_label
    so_accession
    so_name
    feature_class
    feature_name
    file
    file_type
  );
  
  my (
    $dbID,
    $name,
    $logic_name,
    $description,
    $display_label,
    $so_accession,
    $so_name,
    $feature_class,
    $feature_name,
    $file,
    $file_type,
  )
    = rearrange([ @field ], @_);

  $self->dbID            ($dbID);
  $self->name            ($name);
  $self->logic_name      ($logic_name);
  $self->description     ($description);
  $self->display_label   ($display_label);
  $self->so_accession    ($so_accession);
  $self->so_name         ($so_name);
  $self->feature_class   ($feature_class);
  $self->feature_name    ($feature_name);
  $self->file            ($file);
  $self->file_type       ($file_type);

  return $self;
}

sub dbID           { return shift->_generic_get_or_set('dbID',            @_) }
sub name           { return shift->_generic_get_or_set('name',            @_) }
sub logic_name     { return shift->_generic_get_or_set('logic_name',      @_) }
sub description    { return shift->_generic_get_or_set('description',     @_) }
sub display_label  { return shift->_generic_get_or_set('display_label',   @_) }
sub so_accession   { return shift->_generic_get_or_set('so_accession',    @_) }
sub so_name        { return shift->_generic_get_or_set('so_name',         @_) }
sub feature_class  { return shift->_generic_get_or_set('feature_class',   @_) }
sub feature_name   { return shift->_generic_get_or_set('feature_name',    @_) }
sub file           { return shift->_generic_get_or_set('file',            @_) }
sub file_type      { return shift->_generic_get_or_set('file_type',       @_) }

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


