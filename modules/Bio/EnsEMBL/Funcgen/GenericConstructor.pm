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

package Bio::EnsEMBL::Funcgen::GenericConstructor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );

use Role::Tiny;

sub new {
  my ($class, @args) = @_;
  my $self = bless {}, $class;

  $self->_initialise_fields(@args);
  $self->init(@args);
  return $self;
}

# This is a hook for the the object consuming this role to do more
# specific kinds of initialisation.
#
sub init {}

# This is a hook for the the object consuming this role to do define
# the constructor parameters and the setters they should be forwarded
# to.
#
sub _constructor_parameters {
  return {};
}

sub _initialise_fields {

  my $self = shift;
  my @parameters = @_;
  
  my $constructor_key_to_set_method = $self->_constructor_parameters;
  
  if (ref $constructor_key_to_set_method ne 'HASH') {
    throw("_constructor_parameters in " . (ref $self) . " must return a hash reference!");
  }
  
  my @accepted_constructor_parameters = keys %$constructor_key_to_set_method;
  
  return if !@accepted_constructor_parameters;
#   use Data::Dumper;
#   print Dumper(\@accepted_constructor_parameters);

  my @value = rearrange([ @accepted_constructor_parameters ], @parameters);
  
  for (my $index = 0; $index<@value; $index++) {
  
    my $constructor_parameter_name = $accepted_constructor_parameters[$index];
    my $value_to_set               = $value[$index];
    
    my $setter_method = $constructor_key_to_set_method->{$constructor_parameter_name};
    
    if (defined $value_to_set) {
      
      if (! $self->can($setter_method)) {
        throw(
            "The setter method \"" . $setter_method . "\""
            . " to store the parameter \"" . $constructor_parameter_name . "\""
            . " in the constructor of " . (ref $self) 
            . " has not been implemented!");
      }
      $self->$setter_method($value_to_set);
    }
  }
  return;
}

1;
