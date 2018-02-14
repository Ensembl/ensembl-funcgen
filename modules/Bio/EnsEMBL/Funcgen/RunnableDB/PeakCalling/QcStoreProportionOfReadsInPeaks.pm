package Bio::EnsEMBL::Funcgen::RunnableDB::PeakCalling::QcStoreProportionOfReadsInPeaks;

use warnings;
use strict;
use Data::Dumper;

use base 'Bio::EnsEMBL::Hive::Process';

sub run {
  my $self = shift;
  
  my $peak_calling_name            = $self->param_required('peak_calling');
  my $species                      = $self->param_required('species');
  my $num_reads_in_peaks           = $self->param_required('num_reads_in_peaks');
  my $num_reads_in_total           = $self->param_required('num_reads_in_total');
  my $proportion_of_reads_in_peaks = $self->param_required('proportion_of_reads_in_peaks');

  my $peak_calling_adaptor = Bio::EnsEMBL::Registry
    ->get_adaptor(
        $species, 
        'funcgen', 
        'PeakCalling'
    );
  my $peak_calling = $peak_calling_adaptor->fetch_by_name($peak_calling_name);

  my $frip_adaptor = Bio::EnsEMBL::Registry
    ->get_adaptor(
        $species, 
        'funcgen', 
        'frip'
    );
  
  use Bio::EnsEMBL::Funcgen::Frip;
  my $frip = Bio::EnsEMBL::Funcgen::Frip->new(
    -frip            => $proportion_of_reads_in_peaks,
    -total_reads     => $num_reads_in_total,
    -peak_calling_id => $peak_calling->dbID
  );
  $frip_adaptor->store($frip);
  return;
}

1;
