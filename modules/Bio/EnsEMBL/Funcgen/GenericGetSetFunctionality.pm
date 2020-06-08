=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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
use Data::Dumper;

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
  if (! $self->can($dbID_method)) {
    throw(ref($self) . " is missing the method $dbID_method");
  }
  
  my $object_id = $self->$dbID_method;
  if (! defined $object_id) {
    use Carp;
    confess("The method $dbID_method returned an undefined value!");
  }
  my $object_adaptor;
  
  eval {
    $object_adaptor = $self->adaptor;
  };
  if ($@) {
    my ($package, $filename, $line) = caller;
    warn("Calling adaptor failed!");
    warn("\tTrying the older 'db' method.");
    $object_adaptor = $self->adaptor;
    if ($object_adaptor) {
      warn("\tThe db method succeeded.\n");
      warn("\tPlease switch to using 'adaptor' instead of 'db' in $package ($filename).\n");
    }
  }
  
  if (! defined $object_adaptor) {
      throw(
        "Can't get db adaptor from\n" 
        . Dumper($self) . "\n"
        . "This can happen, if " . ref($self) . " does not have db and dbID as constructor parameters.\n"
        . "The method should include them:\n"
        . "\n"
        . "\t" . "sub _constructor_parameters {\n"
        . "\t" . "  return {\n"
        . "\t" . "    dbID           => 'dbID',\n"
        . "\t" . "    db             => 'db',\n"
        . "\t" . "    ... # everything else here\n"
        . "\t" . "  }\n"
      );
  }
  
  my $dba = $object_adaptor->db;
  my $specific_adaptor = $dba->$get_adaptor_method_name;
  my $object = $specific_adaptor->fetch_by_dbID($object_id);
  $self->{$hash_key} = $object;
  return $object;
}

sub _generic_get_or_set {
  my $self  = shift;
  my $name  = shift;
  my $value = shift;
  my $unset = shift;

  if(defined $value) {
    $self->{$name}  = $value;
  }
  if($unset) {
    $self->{$name}  = undef;
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

  if (defined $type and ! $obj->isa($type)) {
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


