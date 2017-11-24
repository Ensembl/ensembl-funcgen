package Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::Constants;

=head1 Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::Constants

Here is where all the constants are declared. To use them, add:

  use Bio::EnsEMBL::Funcgen::ChIPSeqAnalysis::Constants qw ( :all );
  
to your module, then use:
  
  if ($peak_calling_strategy eq CALL_BROAD_PEAKS ) {
    convert_bam_to_bed();
  }

=cut

use base qw( Exporter );
use vars qw( @EXPORT_OK );
use strict;

our @EXPORT_OK = qw(

  CALL_BROAD_PEAKS  
  CALL_NARROW_PEAKS
  
  SKIP_IDR
  RUN_IDR_ON_BIOLOGICAL_REPLICATES 
  RUN_IDR_ON_TECHNICAL_REPLICATES  
  
  BAM_FORMAT
  BIGWIG_FORMAT
  BED_FORMAT
);

our %EXPORT_TAGS = (all => \@EXPORT_OK);

use constant {

  CALL_BROAD_PEAKS  => 'CALL_BROAD_PEAKS',
  CALL_NARROW_PEAKS => 'CALL_NARROW_PEAKS',

  SKIP_IDR                         => 'SKIP_IDR',
  RUN_IDR_ON_BIOLOGICAL_REPLICATES => 'RUN_IDR_ON_BIOLOGICAL_REPLICATES',
  RUN_IDR_ON_TECHNICAL_REPLICATES  => 'RUN_IDR_ON_TECHNICAL_REPLICATES',

  BAM_FORMAT    => 'bam',
  BIGWIG_FORMAT => 'bigwig',
  BED_FORMAT    => 'bed',
};

1;
