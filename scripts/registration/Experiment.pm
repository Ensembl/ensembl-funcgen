package Experiment;
use strict;
use warnings;
use Moose;
use Moose::Util::TypeConstraints;

my $meta = __PACKAGE__->meta;

has 'accession' => (
  is  => 'ro',
  isa => 'Str',
  required => 1,
);

has 'assay_term_id' => (
  is  => 'ro',
  isa => 'Str',
  required => 1,
);

has 'assay_term_name' => (
  is  => 'ro',
  isa => 'Str',
  required => 1,
);

has 'feature_type' => (
  is  => 'ro',
  isa => 'Str',
  );

has 'Experiment' => (
  is  => 'rw',
  isa => 'ArrayRef[Experiment]'
);

has 'File' => (
  is => 'rw',
  isa => 'ArrayRef[File]'
);

1;