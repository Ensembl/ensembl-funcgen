package Bio::EnsEMBL::Funcgen::Hive::Ftp::SplitString;

use strict;
use Data::Dumper;
use base ('Bio::EnsEMBL::Hive::Process');

sub run {
  my $self      = shift;
  my $string         = $self->param('string');
  my $separator      = $self->param('separator');
  my $list_item_name = $self->param('list_item_name');

  my @list_item = split $separator, $string;
  
  foreach my $current_list_item (@list_item) {
    $self->dataflow_output_id({ $list_item_name => $current_list_item }, 2);
  }
}

1;
