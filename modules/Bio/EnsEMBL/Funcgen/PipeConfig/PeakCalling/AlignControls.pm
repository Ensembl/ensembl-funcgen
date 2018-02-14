package Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::AlignControls;

use strict;
use warnings;
use base 'Bio::EnsEMBL::Funcgen::PipeConfig::PeakCalling::Base';
use Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf;

sub pipeline_analyses {
    my $self = shift;
    return [
        @{
            $self->generate_parallel_alignment_analyses({
                start  => 'start_align_controls',
                prefix => '',
                suffix => '_controls',
            })
        },
    ];
}

1;

