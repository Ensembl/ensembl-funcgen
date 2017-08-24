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
use Bio::EnsEMBL::Utils::Scalar    qw( assert_ref );
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );

=head2 new
  Example      : 
  Description : Constructor method
  Returntype  : Bio::EnsEMBL::Funcgen::ReadFileExperimentalConfiguration
  Exceptions  : None 
  Caller      : General
  Status      : Stable

=cut

sub new {
  my $caller = shift;
  my $class  = ref($caller) || $caller;
  my $self   = bless {}, $class;
  
  my @parameters = @_;
  
  my @simple_accessor_fields = $self->_simple_accessor_fields;
  my @setter_fields          = $self->_setter_fields;
  
  push @simple_accessor_fields, 'db';
  push @simple_accessor_fields, 'dbID';
  
  $self->init_fields(\@parameters, \@simple_accessor_fields);
  $self->init_fields(\@parameters, \@setter_fields, 'set_');

  return $self;
}

sub dbID                 { return shift->_generic_get_or_set('dbID', @_) }
sub db                   { return shift->_generic_get_or_set('db',   @_) }

sub init_fields {

  my $self          = shift;
  my $parameter     = shift;
  my $init_fields   = shift;
  my $method_prefix = shift;
  
  if (! defined $method_prefix) {
    $method_prefix = '';
  }
  
  my @value = rearrange([ @$init_fields ], @$parameter);
  
  for (my $index = 0; $index<@value; $index++) {
  
    my $setter_method;
    
    if ($method_prefix) {
      my $field_name = $init_fields->[$index];
      
      # Remove underscores and replace them by camel case convention to 
      # get the object names.
      #
      $field_name =~ s/_(.)/uc($1)/ge;
      $field_name =~ s/(^.)/uc($1)/ge;
      $setter_method = $method_prefix . $field_name;
    } else {
      $setter_method = $init_fields->[$index];
    }
    
    if (! $self->can($setter_method)) {
      throw("Unknown method " . $setter_method . "!");
    }
    my $value_to_set  = $value[$index];
    if (defined $value_to_set) {
      $self->$setter_method($value_to_set);
    }
  }
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
  assert_ref($obj, $type);
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


