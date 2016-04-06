=head1 NAME

Bio::EnsEMBL::Funcgen::Hive::ErsaCleanup

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::Hive::ErsaCleanup;

use base ('Bio::EnsEMBL::Funcgen::Hive::BaseDB');

use warnings;
use strict;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd );

sub run {
  my $self = shift; 
  
  # Removes .sav files from the peaks directory. These are saved R sessions
  # for which there is no use case.
  #
  # Looks like RunIDR writes these and PostProcessIDRReplicates reads them, 
  # so we have to be selective and use strain names here to narrow this down
  # to the .sav files belonging to cell lines that are actually done.
  #
  # my $peaks_output_dir = $self->peaks_output_dir;
  # my $rm_sav_cmd = qq(bash -o pipefail -c "find $peaks_output_dir -name '*.sav' | xargs rm -f ");
  # run_system_cmd($rm_sav_cmd);
  
  my $file_to_delete = $self->param('file_to_delete');
  
  foreach my $current_file (@$file_to_delete) {
    my $cmd = "rm -f $current_file";
    print $cmd . "\n";
    run_system_cmd($cmd); 
  }
  return;
}

1;
