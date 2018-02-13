package Bio::EnsEMBL::Funcgen::Utils::Fastq::Parser;

use strict;
use Data::Dumper;
use Carp;

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

sub parse_next_record {

  my $self = shift;
  my $fh   = shift;
  
  if (! defined $fh) {
    confess("File handle is not defined!");
  }

  my @current_fastq_record;
  
  foreach my $line_of_record (1, 2, 3, 4) {
    my $line = $fh->getline;
    
#     print "line = $line\n";
#     print Dumper(\@current_fastq_record);
    
    if (! defined $line) {
      return;
    }
    
    chomp $line;
    push @current_fastq_record, $line;
  }
  return \@current_fastq_record;
}

1;
