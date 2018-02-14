package Bio::EnsEMBL::Funcgen::PeakCallingPlan::select_EnsemblAlignmentAnalysis;

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw ( :all );
use Carp;
use Role::Tiny;

sub select_EnsemblAlignmentAnalysis {

  my $param = shift;

  my $experiment = $param->{experiment};
  my $has_single_end_reads_only = $param->{has_single_end_reads_only};
  my $has_paired_end_reads_only = $param->{has_paired_end_reads_only};

  my $has_mix_of_read_types = ! $has_single_end_reads_only && ! $has_paired_end_reads_only;

  if ($has_mix_of_read_types) {
    #confess("Experiment ". $experiment->name ." has a mix of read types!");
    print STDERR ("Experiment ". $experiment->name ." has a mix of read types!");
    return ENSEMBL_HODGEPODGE_ALIGNMENT_ANALYSIS;
  }
  
  if ($has_single_end_reads_only) {
    return ENSEMBL_SINGLE_END_ALIGNMENT_ANALYSIS;
  }
  if ($has_paired_end_reads_only) {
    return ENSEMBL_PAIRED_END_ALIGNMENT_ANALYSIS;
  }
  
  confess("has_single_end_reads_only and has_paired_end_reads_only can't both be true!");
}

1;
