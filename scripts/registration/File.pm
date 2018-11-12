package File;
use strict;
use warnings;
use Moose;
use Moose::Util::TypeConstraints;
 

has 'accession' => (
  is => 'ro',
  isa => 'Str',
  required => 1,
    
);

has 'file_type' => (
  is        => 'ro',
  isa       => enum( [ qw[ fastq ]] ),
  required  => 1,
);

# has 'submitted_file_name' => (
#   is => 'ro',
#   isa => 'Str',
#   required => 1,
# );

has 'md5sum' => (
  is        => 'ro',
  isa       => 'Str',
  required  => 0,
);

has 'read_length' => (
  is        => 'ro',
  isa       => 'Int',
  required  => 0,
);

has 'read_count' => (
  is        => 'ro',
  isa       => 'Int',
  required  => 0,
);

has 'run_type' => (
  is => 'ro',
  isa => enum([qw[ single-ended paired-ended ]]),
  required => 1,
);

has 'aliases' => (
  is      => 'rw',
  isa     => 'ArrayRef[Str]',
  traits  => [ 'Array' ],
  handles => {
    file_push     => 'push',
    file_shift    => 'shift',
    file_elements => 'elements',
    file_count    => 'count',
    file_is_empty => 'is_empty',
    },
  required => 0,
);

# has 'gender' => (
#   is => 'ro',
#   isa => enum([qw[ male female hermaphrodite mixed unknown ]]),
#   required => 1
# );


1;
