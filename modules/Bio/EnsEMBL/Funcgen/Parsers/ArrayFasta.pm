package Bio::EnsEMBL::Funcgen::Parsers::ArrayFasta;

use strict;
use Data::Dumper;

sub new {
  my ( $class, @args ) = @_;
  my $self = bless {}, $class;  
  return $self;
}

sub validate_configuration {
  my $self  = shift;
  my $param = shift;
  
  my $IIDREGEXP   = $param->{IIDREGEXP};
  my $IFIELDORDER = $param->{IFIELDORDER};
  my $probe_file  = $param->{probe_file};

  my %valid_fields = (
    -probe_set   => undef,
    -name        => undef,
    -array       => undef,
    -array_chip  => undef,
    -description => undef,
  );
  my %field_order = %$IFIELDORDER;

  foreach my $config_field (keys(%field_order)) {

    if (! exists $valid_fields{$config_field}) {
      throw("Found invalid field on ImportArrays.pm config:\t$config_field\n".
            "IFIELDORDER must only contain keys:\t".join("\t", keys %valid_fields));
    }
  }
}

sub run {

  my $self  = shift;
  my $param = shift;
  
  my $header_regex = $param->{IIDREGEXP};
  my $field_order  = $param->{IFIELDORDER};
  my $probe_file   = $param->{probe_file};
  my $call_back    = $param->{call_back};
 
  # This is done so we can dynamically assign fields using their match order
  # Don't actually use 4 and 5 here, but here just in case config regexps are 
  # extended.
  #
  my @match_refs = (\$1, \$2, \$3, \$4, \$5);
  
  use Bio::EnsEMBL::IO::Parser::Fasta;
  my $fasta_parser = Bio::EnsEMBL::IO::Parser::Fasta->open($probe_file);
  
  while ($fasta_parser->next) {

    my $header   = $fasta_parser->getHeader;
    my $sequence = $fasta_parser->getSequence;
    
    my $header_as_our_regex_expects_it = '>' . $header;

    # This places header values into @match_refs
    if ($header_as_our_regex_expects_it !~ /$header_regex/) {
      use Carp;
      confess("Problem with header: ${header_as_our_regex_expects_it}, it does not match the regular expression ${header_regex}!");
    }
    
    if ($sequence =~ /^>/) {
      confess("Problem with sequence ($sequence), it starts with \">\", so probably a problem parsing the fasta file!");
    }
    
    my $probe_data = {
      -sequence => $sequence
    };
    
    foreach my $field (keys %$field_order) {
      $probe_data->{$field} = ${$match_refs[$field_order->{$field}]};
    }
    $call_back->($probe_data);
  }
  return;
}

1;

