package Bio::EnsEMBL::Funcgen::Parsers::DataDumper;

use strict;
use Data::Dumper;

my $record_terminator = ";\n";

sub new {
  my ( $class, @args ) = @_;
  my $self = bless {}, $class;  
  return $self;
}

=head2 load_first_item_from_data_dump_file

  Convenience method, if there is only one item and it is small.

=cut
sub load_first_item_from_data_dump_file {
  my $self = shift;
  my $file = shift;
  my $array_ref = $self->load_data_dump_file($file);
  return $array_ref->[0];
}

=head2 load_data_dump_file

  Convenience method, if the data is small.

=cut
sub load_data_dump_file {
  my $self = shift;
  my $file = shift;
  my @result;

  $self->parse({
    data_dumper_file => $file,
    call_back        => sub {
      push @result, shift;
    },
  });
  return \@result;
}

sub parse {

  my $self  = shift;
  my $param = shift;
  
  my $data_dumper_file = $param->{data_dumper_file};
  my $call_back        = $param->{call_back};
  
  if (! -e $data_dumper_file) {
    use Bio::EnsEMBL::Utils::Exception qw( throw );
    throw("The file $data_dumper_file does not exist!");
  }
  local $/ = $record_terminator;
  my $in;
  open $in, $data_dumper_file;

  LINE:
  while (my $current_record = <$in>) {

    # Skip whitespace only records. This should only occurr at the end of the 
    # file.
    next LINE if ($current_record =~ /^\s+$/m);

    no strict;
    my $parsed = eval $current_record;
    use strict;
    
    eval {
      $call_back->($parsed);
    };
    if ($@) {
      use Carp;
      confess(
        "Error processing " . $current_record . "!\n"
        . "The error message is: " . $@
      );
    }
  }
  $in->close;
}

1;

