package Bio::EnsEMBL::Funcgen::Config::ArrayInstance;

use strict;
use Data::Dumper;
use List::Util qw(any);
use Bio::EnsEMBL::Funcgen::Config::ImportArrays;
use Bio::EnsEMBL::Analysis::Tools::Utilities qw(parse_config);

sub new {
  my ( $class, @args ) = @_;
  my $self = bless {}, $class;
  return $self;
}

sub read_and_check_config {

  my $self = shift;
  my $array_format = shift;

  my $former_logic_name = 'IMPORT_'. uc($array_format) .'_ARRAYS';
  parse_config($self, $ARRAY_CONFIG, $former_logic_name);

  foreach my $config_var (qw( IIDREGEXP IFIELDORDER INPUT_FORMAT )) {
    if ( ! defined $self->$config_var ) {
      throw("You must define $config_var.");
    }
  }
}

sub get_ARRAY_PARAMS_by_array_name {
  my ( $self, $array_name ) = @_;

  if (any { $_ eq $array_name } @{$self->{'_CONFIG_ARRAYS_WITH_DEFAULT_PARAMS'}})
  {
    $self->{'_CONFIG_ARRAY_PARAMS'}->{$array_name} = $self->{'_CONFIG_ARRAY_PARAMS'}->{'Default'};
    $self->{'_CONFIG_ARRAY_PARAMS'}->{$array_name}->{'-name'} = $array_name;
   } elsif(! exists $self->{'_CONFIG_ARRAY_PARAMS'}{$array_name}) {
    use Carp;
    confess("No ARRAY_PARAMS config available for $array_name.  You must add this to the ImportArrays config before importing");
 }

  return $self->{'_CONFIG_ARRAY_PARAMS'}{$array_name};
}

sub IIDREGEXP {
  my ( $self, $value ) = @_;
  $self->{'_CONFIG_IIDREGEXP'} = $value  if defined $value;
  return $self->{'_CONFIG_IIDREGEXP'};
}

sub get_IIDREGEXP{
  my $self = shift;
  return $self->{'_CONFIG_IIDREGEXP'};
}

sub get_IFIELDORDER{
  my $self = shift;
  return $self->{'_CONFIG_IFIELDORDER'};
}

sub ARRAY_FORMAT {
  my ( $self, $value ) = @_;
  $self->{'_ARRAY_FORMAT'} = $value if defined $value;
  return $self->{'_ARRAY_FORMAT'};
}

sub INPUT_FORMAT{
  my ( $self, $value ) = @_;
  $self->{'_CONFIG_INPUT_FORMAT'} = $value  if defined $value ;
  return $self->{'_CONFIG_INPUT_FORMAT'};
}

sub ARRAY_PARAMS {
  my ( $self, $value ) = @_;
  $self->{'_CONFIG_ARRAY_PARAMS'} = $value  if defined $value;
  return $self->{'_CONFIG_ARRAY_PARAMS'};
}

sub ARRAYS_WITH_DEFAULT_PARAMS {
  my ( $self, $value ) = @_;
  $self->{'_CONFIG_ARRAYS_WITH_DEFAULT_PARAMS'} = $value  if defined $value;
  return $self->{'_CONFIG_ARRAY_PARAMS'};
}

sub IFIELDORDER {
  my ( $self, $value ) = @_;
  $self->{'_CONFIG_IFIELDORDER'} = $value if defined $value;
  return $self->{'_CONFIG_IFIELDORDER'};
}

1;
