=head1 NAME

Bio::EnsEMBL::Funcgen::Hive::ErsaCleanup

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::Hive::ErsaCleanup;

use base ('Bio::EnsEMBL::Funcgen::Hive::Base');

use warnings;
use strict;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd );

sub run {
  my $self = shift; 
  
  my $file_to_delete_after_cell_line_has_been_processed = $self->param('file_to_delete_after_cell_line_has_been_processed');
  
  foreach my $current_file (@$file_to_delete_after_cell_line_has_been_processed) {
    my $cmd = "rm -f $current_file";
    print $cmd . "\n";
    run_system_cmd($cmd); 
  }
  return;
}

1;
