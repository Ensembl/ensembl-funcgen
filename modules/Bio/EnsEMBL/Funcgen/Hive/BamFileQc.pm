package Bio::EnsEMBL::Funcgen::Hive::BamFileQc;

use warnings;
use strict;

use base qw (Bio::EnsEMBL::Hive::Process );

sub run {
  my $self = shift;
   
  my $bam_file = $self->param('bam_file_for_qc');
  
  if (! -e $bam_file) {
  
    use Bio::EnsEMBL::Utils::Logger;
    my $logger = Bio::EnsEMBL::Utils::Logger->new();
    $logger->error(
      "The bam file $bam_file doesn't exist. This can happen, if it has just "
      . "been created and the file system has not been updated yet on all "
      . "nodes. In cases like these it may be resolved by restarting after "
      . "waiting a bit. This job will sleep for 10 seconds now and then die. "
      . "Hopefully upon retry it will work.\n"
    );
    sleep(10);
    die("$bam_file doesn't exist. Dying now and hoping for more luck in "
      . "the next life.");
  }
  
  return;
}

1;
