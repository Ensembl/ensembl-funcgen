package Bio::EnsEMBL::Funcgen::Utils::ERSAShared;

use warnings;
use strict;

use base qw( Exporter );
use vars qw( @EXPORT_OK );

@EXPORT_OK = qw(
  create_replicate_input_subset_string
);

sub create_replicate_input_subset_string {

  my @input_subset_replicates = @_;

  my @sorted_replicate_number = sort {
       $a->biological_replicate <=> $b->biological_replicate
    || $a->technical_replicate  <=> $b->technical_replicate
  } @input_subset_replicates;
  
  my @replicate_number_strings;
  foreach my $current_replicate_number (@sorted_replicate_number) {
    push @replicate_number_strings,
        'BR' . $current_replicate_number->biological_replicate
      . 'TR' . $current_replicate_number->technical_replicate;
  }
  my $replicate_number_string = join('_', sort(@replicate_number_strings));
  
  return $replicate_number_string;
}

1;
