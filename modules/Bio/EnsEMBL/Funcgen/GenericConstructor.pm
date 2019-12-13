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

package Bio::EnsEMBL::Funcgen::GenericConstructor;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Data::Dumper;
#use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw warning deprecate );

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
  
  my $class = ref $self;
  
  if (ref $constructor_key_to_set_method ne 'HASH') {
    throw("_constructor_parameters in " . (ref $self) . " must return a hash reference!");
  }
  
  my @accepted_constructor_parameters = keys %$constructor_key_to_set_method;
  
  return if !@accepted_constructor_parameters;

  my @value = rearrange_pp([ @accepted_constructor_parameters ], @parameters);
  
  my %all_parameters = @parameters;
  
  for (my $index = 0; $index<@value; $index++) {
  
    my $constructor_parameter_name = $accepted_constructor_parameters[$index];
    my $value_to_set               = $value[$index];
    
    my $setter_method = $constructor_key_to_set_method->{$constructor_parameter_name};
    
    if ($setter_method eq '') {
        throw(
            "\n\tError parsing the constructor parameters. This usually"
            . " happens, \n\tif the parameter names in the constructor weren't"
            . " prefixed \n\twith a dash:\n"
            . "Wrong:\n"
            . "\t$class->new( foo => bar )\n"
            . "Correct:\n"
            . "\t$class->new( -foo => bar )\n"
        );
    }
    
    # Tick the parameters that have been used off the list
    delete $all_parameters{'-' . $constructor_parameter_name};

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
  
  # Don't complain about these until all API objects use the new 'adaptor' 
  # name.
  #
  my @unchecked_parameters = (
    '-adaptor',
    '-db'
  );
  
  my @unused_parameters;
  PARAMETER:
  foreach my $current_parameter (keys %all_parameters) {
  
    my $is_unchecked_parameter = grep { $current_parameter ne $_ } @unchecked_parameters;
  
    if ($is_unchecked_parameter) {
      next PARAMETER;
    }
    push @unused_parameters, $current_parameter;
  }
  
  if (@unused_parameters) {
    
    my %as_hash = @parameters;
    $Data::Dumper::Maxdepth = 2;
    throw(
       "The following parameters were passed as parameters to the constructor of\n\n"
       . "\t$class\n\n"
     . "but they are not used:\n\n" 
     . Dumper(\%all_parameters) . "\n\n"
     . "The complete set of parameters passed to the constructor was:\n\n"
     . Dumper( \%as_hash ) . "\n\n"
     . "The constructor accepts the following parameters:\n\n"
     . Dumper( [ map { '-' . $_ } sort @accepted_constructor_parameters ] )
    );
  }
  return;
}

# Using pure perl implementation copied over from 
# Bio::EnsEMBL::Utils::Argument until
# https://www.ebi.ac.uk/panda/jira/browse/ENSCORESW-2464
# has been resolved.
#
sub rearrange_pp {
  my $order = shift;

  if ( $order eq "Bio::EnsEMBL::Utils::Argument" ) {
    # skip object if one provided
    $order = shift;
  }

  # If we've got parameters, we need to check to see whether
  # they are named or simply listed. If they are listed, we
  # can just return them.
  unless ( @_ && $_[0] && substr( $_[0], 0, 1 ) eq '-' ) {
    return @_;
  }
  
  # Push undef onto the end if % 2 != 0 to stop warnings
  push @_,undef unless $#_ %2;
  my %param;
  while( @_ ) {
    #deletes all dashes & uppercases at the same time
    (my $key = shift) =~ tr/a-z\055/A-Z/d;
    $param{$key} = shift;
  }
  
  # What we intend to do is loop through the @{$order} variable,
  # and for each value, we use that as a key into our associative
  # array, pushing the value at that key onto our return array.
  return map { $param{uc($_)} } @$order;
}

1;
