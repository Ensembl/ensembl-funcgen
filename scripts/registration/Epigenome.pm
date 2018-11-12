package Epigenome;
use strict;
use warnings;
use Moose;
use Moose::Util::TypeConstraints;

my $meta = __PACKAGE__->meta;

has 'accession' => (
  is => 'ro',
  isa => 'Str',
  required => 1,
);



has 'biosample_term' => (
  is => 'ro',
  isa => 'Str',
  required => 0,
);

has 'assay_type' => (
  is => 'ro',
  isa => 'Str',
  required => 0,
);

has 'description' => (
  is => 'ro',
  isa => 'Str',
  required => 0,
);

has 'project' => (
  is => 'ro',
  isa => enum([qw[ encode blueprint ddbj roadmap ]]),
  required => 1,
);

has 'status' => (
  is => 'ro',
  isa => 'Str',
  required => 0,
);

has 'species' => (
  is => 'ro',
  isa => 'Str',
  required => 0,
);

has 'ihec' => (
  is => 'ro',
  isa => 'Str',
  required => 0,
);

has 'assay_term_id' => (
  is => 'ro',
  isa => 'Str',
  required => 0,
);

has 'biosample_type' => (
  is => 'ro',
  isa => 'Str',
  required => 0,
);

has 'assay_synonyms' => (
  is => 'ro',
  isa => 'Str',
  required => 0,
);

has 'biosample_term_id' => (
  is => 'ro',
  isa => 'Str',
  required => 0,
);

has 'biosample_synonyms' => (
  is => 'ro',
  isa => 'Str',
  required => 0,
);



# has 'gender' => (
#   is => 'ro',
#   isa => enum([qw[ male female hermaphrodite mixed unknown ]]),
#   required => 1
# );

has 'Experiment' => (
  is => 'rw',
  isa => 'ArrayRef[Experiment]'
);
# https://www.encodeproject.org/reference-epigenomes/ENCSR649KUM/?format=json
#

sub print_all_attributes {
  my ($self) = @_;

  for my $attr ( $meta->get_all_attributes ) {
    print $attr->name .":\t";
    print $self->{$attr->name}, "\n" if(defined $self->{$attr->name});
  }
}





1;