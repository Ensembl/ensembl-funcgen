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

sub new {
  my $caller = shift;
  my $class  = ref($caller) || $caller;
  my $self   = bless {}, $class;
#   die "Test!\n";
  $self->init(@_);
  return $self;
}

sub init {
  my $self = shift;
  $self->_build_class_methods(@_);
}

sub _build_class_methods {

  my $self = shift;
  my @parameters = @_;
  
  my $constructor_key_to_set_method = $self->_constructor_parameters;
  
  if (ref $constructor_key_to_set_method ne 'HASH') {
    throw("_constructor_parameters in " . (ref $self) . " must return a hash reference!");
  }
  
  $constructor_key_to_set_method->{db}   = 'db';
  $constructor_key_to_set_method->{dbID} = 'dbID';
  
  my $package = ref $self;
  my $is_constructed_flag = '__is_constructed';
  
  if (! defined $self->can($is_constructed_flag)) {
  
    $self->_build_simple_accessors;
    $self->_build_fetch_methods;
    $self->_build_get_methods;
    $self->_build_set_methods;
    
    no strict;
    *{$package . '::' . $is_constructed_flag} = sub { return 1 };
    use strict;
  }
  $self->_init_fields(\@parameters, $constructor_key_to_set_method);
}

sub _constructor_parameters { return {} }
sub _simple_accessors       { return [] }
sub _fetch_methods          { return [] }
sub _get_methods            { return [] }
sub _set_methods            { return [] }

sub dbID { return shift->_generic_get_or_set('dbID', @_) }
sub db   { return shift->_generic_get_or_set('db',   @_) }

sub _init_fields {

  my $self          = shift;
  my $parameter     = shift;
  my $constructor_key_to_set_method = shift;
  
  my @accepted_constructor_parameters = keys %$constructor_key_to_set_method;

  my @value = rearrange([ @accepted_constructor_parameters ], @$parameter);
  
  for (my $index = 0; $index<@value; $index++) {
  
    my $constructor_parameter_name = $accepted_constructor_parameters[$index];
    my $value_to_set               = $value[$index];
    
    my $setter_method = $constructor_key_to_set_method->{$constructor_parameter_name};
    
    if (defined $value_to_set) {
      
      if (! $self->can($setter_method)) {
        throw("The setter method \"" . $setter_method . "\" to store the parameter \"" . $constructor_parameter_name . "\" in the constructor of " . (ref $self) . " has not been implemented!");
      }
      $self->$setter_method($value_to_set);
    }
  }
  return;
}

sub _build_simple_accessors {
  my $self = shift;
  my $fetch_method_specification = $self->_simple_accessors;
  
  foreach my $current_fetch_method_specification (@$fetch_method_specification) {
    $self->_build_simple_accessor($current_fetch_method_specification)
  }
  return;
}

sub _build_simple_accessor {
  my $self  = shift;
  my $specs = shift;
  
  my $method_name = $specs->{method_name};
  my $hash_key    = $specs->{hash_key};

  my $package = ref $self;
  my $full_method_name = $package . "::" . $method_name;
  
  no strict;
  
  *{$full_method_name} = sub { return shift->_generic_get_or_set($hash_key, @_) };
  use strict;

  return;
}

sub _build_get_methods {
  my $self = shift;
  my $get_method_specification = $self->_get_methods;
  
  foreach my $current_get_method_specification (@$get_method_specification) {
    $self->_build_get_method($current_get_method_specification)
  }
  return;
}

sub _build_get_method {
  my $self  = shift;
  my $specs = shift;
  
  my $method_name = $specs->{method_name};
  my $hash_key    = $specs->{hash_key};

  my $package = ref $self;
  
  no strict;
  *{$package . "::" . $method_name} = sub { return shift->_generic_get($hash_key,  @_) };
  use strict;

  return;
}

sub _build_set_methods {
  my $self = shift;
  my $set_method_specification = $self->_set_methods;
  
  foreach my $current_set_method_specification (@$set_method_specification) {
    $self->_build_set_method($current_set_method_specification)
  }
  return;
}

sub _build_set_method {
  my $self  = shift;
  my $specs = shift;
  
  my $method_name   = $specs->{method_name};
  my $hash_key      = $specs->{hash_key};
  my $expected_type = $specs->{expected_type};

  my $package = ref $self;
  
  no strict;
  *{$package . "::" . $method_name} = sub { 
    return shift->_generic_set(
      $method_name, $hash_key, $expected_type, @_
    )
  };

  use strict;

  return;
}

sub _build_fetch_methods {
  my $self = shift;
  my $fetch_method_specification = $self->_fetch_methods;  
  foreach my $current_fetch_method_specification (@$fetch_method_specification) {
    $self->_build_fetch_method($current_fetch_method_specification)
  }
  return;
}
sub _build_fetch_method {

  my $self  = shift;
  my $specs = shift;
  
  
  my $method_name             = $specs->{method_name};
  my $hash_key                = $specs->{hash_key};
  my $get_adaptor_method_name = $specs->{get_adaptor_method_name};
  my $dbID_method             = $specs->{dbID_method};

  my $package = ref $self;  

  no strict;
  
  *{$package . "::" . $method_name} = 
  sub {
      my $self = shift;

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
  };

  use strict;
  return;
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
  my $method_name = shift;
  my $name  = shift;
  my $type  = shift;
  my $obj   = shift;
  
  if (! defined $obj) {
    throw("$name was not defined!");
  }

  if (! $obj->isa($type)) {
    throw("Expected $type, but got " . (ref $obj) . " when calling " . $method_name);
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


