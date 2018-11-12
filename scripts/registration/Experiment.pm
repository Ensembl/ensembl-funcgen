package Experiment;
use strict;
use warnings;
use Moose;
use Moose::Util::TypeConstraints;

extends 'Epigenome';


has 'assay_term_id' => (
  is => 'ro',
  isa => 'Str',
  required => 1,
);

has 'assay_title' => (
  is => 'ro',
  isa => 'Str',
  required => 1,
);

has 'biosample_term_id' => (
  is => 'ro',
  isa => 'Str',
  required => 1,
);

has 'biosample_term_name' => (
  is => 'ro',
  isa => 'Str',
  required => 1,
);

has 'biosample_summary' => (
  is => 'ro',
  isa => 'Str',
  required => 1,
);

has 'biosample_type' => (
  is => 'ro',
  isa => 'Str',
  required => 1,
);

has 'description' => (
  is => 'ro',
  isa => 'Str',
);

1;