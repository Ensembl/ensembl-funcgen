package Bio::EnsEMBL::Funcgen::Parsers::DataDumper;

use strict;
use Data::Dumper;

sub new {
  my ( $class, @args ) = @_;
  my $self = bless {}, $class;  
  return $self;
}

sub parse {

  my $self  = shift;
  my $param = shift;
  
  my $data_dumper_file = $param->{data_dumper_file};
  my $call_back        = $param->{call_back};

  local $/ = ";\n";
  open my $in, $data_dumper_file;

  LINE:
  while (my $current_record = <$in>) {

    # Skip whitespace only records. This should only occurr at the end of the 
    # file.
    next LINE if ($current_record =~ /^\s+$/m);

    no strict;
    my $parsed = eval $current_record;
    use strict;
    
    $call_back->($parsed);
  }
  $in->close;
}

1;

