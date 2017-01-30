package Bio::EnsEMBL::Funcgen::Config::ArrayFormatConfig;

use strict;

sub new {
  my ( $class, @args ) = @_;
  my $self = bless {
    _array_format_config => &_constant_array_format_config(),
  }, $class;
  return $self;
}

sub for_array_class {
  my $self        = shift;
  my $array_class = shift;
  
  if (! defined $array_class) {
    die;
  }

  if (! exists $self->array_format_config->{$array_class}) {
    die("Unknown array class ${array_class}. Please make sure the array has been registered correctly in the Bio::EnsEMBL::Funcgen::Config::ImportArrays module.");
#     warn("Unknown array class ${array_class}. Using default values.");
#     return ArrayClassConfiguration->new($self->default_array_format_config);
  }
  
  return ArrayClassConfiguration->new(
    $self
    ->array_format_config
    ->{$array_class}
  );
}

sub array_format_config {
  my $self        = shift;
  return $self->{_array_format_config};
}

sub default_array_format_config {
  return {
    probeset_arrays      => 0,
    linked_arrays        => 0,
    sense_interrogation  => 0,
  };
}

sub _constant_array_format_config {

  #Array format config for current vendors
  my %array_format_config = (
    AFFY_UTR => {
      probeset_arrays      => 1,
      linked_arrays        => 0,
      sense_interrogation  => 0,
    },
    AFFY_ST => {
      probeset_arrays      => 1,
      linked_arrays        => 1,
      sense_interrogation  => 1,
    },
    ILLUMINA_WG => {
      probeset_arrays      => 0,
      linked_arrays        => 0,
      sense_interrogation  => 0,
    },
    ILLUMINA_INFINIUM => {
      probeset_arrays      => 0,
      linked_arrays        => 0,
      sense_interrogation  => 0,
    },
    AGILENT => {
      probeset_arrays      => 0,
      linked_arrays        => 1,
      sense_interrogation  => 0,
    },
    PHALANX => {
      probeset_arrays      => 0,
      linked_arrays        => 1,
      sense_interrogation  => 0,
    },
    CODELINK => {
      probeset_arrays      => 0,
      linked_arrays        => 1,
      sense_interrogation  => 0,
    },
    LEIDEN => {
      probeset_arrays      => 0,
      linked_arrays        => 1,
      sense_interrogation  => 0,
    },
    STEMPLE_LAB_SANGER => {
      probeset_arrays      => 0,
      linked_arrays        => 1,
      sense_interrogation  => 0,
    },
    NIMBLEGEN_MODENCODE => {
      probeset_arrays      => 0,
      linked_arrays        => 1,
      sense_interrogation  => 0,
    },
    SLRI => {
      probeset_arrays      => 0,
      linked_arrays        => 1,
      sense_interrogation  => 0,
    },
    UCSF => {
      probeset_arrays      => 0,
      linked_arrays        => 1,
      sense_interrogation  => 0,
    },
    WUSTL =>  {
      probeset_arrays      => 0,
      linked_arrays        => 1,
      sense_interrogation  => 0,
    },
  );
  use Hash::Util qw( lock_hash );
  lock_hash(%array_format_config);
  return \%array_format_config;
}

1;

package ArrayClassConfiguration;

use strict;

sub new {
  my ( $class, $array_class_configuration_hash ) = @_;
  my $self = bless $array_class_configuration_hash, $class;
  return $self;
}

sub has_probesets {
  return shift->{probeset_arrays};
}

sub has_linked_arrays {
  return shift->{linked_arrays};
}

sub has_sense_interrogation {
  return shift->{sense_interrogation};
}


1;

