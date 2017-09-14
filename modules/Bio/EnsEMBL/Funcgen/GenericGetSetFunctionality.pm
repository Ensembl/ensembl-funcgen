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

=head1 SYNOPSIS

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw );
use base qw( Exporter );
use vars qw( @EXPORT_OK );

@EXPORT_OK = qw(
  _generic_get_or_set
  _generic_set
  _generic_get
  _generic_fetch
);

sub _generic_fetch {

  my $self                    = shift;
  my $hash_key                = shift;
  my $get_adaptor_method_name = shift;
  my $dbID_method             = shift;

  if ($self->{$hash_key}) {
    return $self->{$hash_key};
  }
  my $object_id = $self->$dbID_method;
    if (! defined $object_id) {
    die;
  }
  my $object_adaptor = $self->db->db->$get_adaptor_method_name;
  my $object = $object_adaptor->fetch_by_dbID($object_id);
  $self->{$hash_key} = $object;
  return $object;
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

sub _generic_set {
  my $self  = shift;
  my $name  = shift;
  my $type  = shift;
  my $obj   = shift;
  
  if (! defined $obj) {
    throw("$name was not defined!");
  }

  if (! $obj->isa($type)) {
    throw("Expected $type, but got " . (ref $obj) . "!");
  }
  
  $self->{$name} = $obj;
  return $obj;
}

sub _generic_get {
  my $self  = shift;
  my $name  = shift;
  
  if ($self->{$name}) {
    return $self->{$name};
  }
  return
}

1;


