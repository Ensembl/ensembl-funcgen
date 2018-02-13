package Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::summarise_ReadFile;

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::Constants qw ( :all );

use Role::Tiny;

sub summarise_ReadFile {

  my $read_file = shift;

  if ($read_file->is_paired_end) {
    my $read_file_mate = $read_file->fetch_mate_ReadFile;
    my $summary = {
        $read_file->paired_end_tag      => $read_file->name,
        $read_file_mate->paired_end_tag => $read_file_mate->name,
        type                            => PAIRED_END,
    };
    return $summary;
  }
  
  my $summary = {
      name => $read_file->name,
      type => SINGLE_END,
  };
  
  return $summary;
}

1;
