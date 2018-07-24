package Bio::EnsEMBL::Funcgen::PeakCallingPlan::PeakCallingPlanBuilder;

use strict;
use Data::Dumper;
use Role::Tiny::With;
use Bio::EnsEMBL::Funcgen::PeakCallingPlan::Constants qw( :all );
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_set
  _generic_get
);

# Input

sub set_signal_alignment        { return shift->_generic_set('signal_alignment',        undef, @_); }
sub set_control_alignment       { return shift->_generic_set('control_alignment',       undef, @_); }
sub set_idr                     { return shift->_generic_set('output_real',             undef, @_); }
sub set_experiment              { return shift->_generic_set('experiment',              undef, @_); }
sub set_peaks_output_dir        { return shift->_generic_set('output_dir',              undef, @_); }
sub set_alignment_namer         { return shift->_generic_set('alignment_namer',         undef, @_); }
sub set_control_alignment_namer { return shift->_generic_set('control_alignment_namer', undef, @_); }
sub set_samtools_fasta_index    { return shift->_generic_set('samtools_fasta_index',    undef, @_); }
sub set_chromosome_lengths_by_species_assembly { 
  return shift->_generic_set(
    'chromosome_lengths_by_species_assembly', 
    undef, 
    @_
  ); 
}

# Output

sub get_peakcalling_plan { 
  return shift->_generic_get('peakcalling');
}

sub get_bam_to_bed_conversion_plan { 
  return shift->_generic_get('bam_to_bed_conversion_plan');
}

sub _set_peakcalling_plan {
  my $self = shift;
  my $obj  = shift;
  return $self->_generic_set('peakcalling', undef, $obj);
}
sub _set_bam_to_bed_conversion_plan {
  my $self = shift;
  my $obj  = shift;
  return $self->_generic_set('bam_to_bed_conversion_plan', undef, $obj);
}


sub _get_signal_alignment        { return shift->_generic_get('signal_alignment');          }
sub _get_control_alignment       { return shift->_generic_get('control_alignment');         }
sub _get_idr                     { return shift->_generic_get('output_real');               }
sub _get_experiment              { return shift->_generic_get('experiment');                }
sub _get_peaks_output_dir        { return shift->_generic_get('output_dir');                }
sub _get_alignment_namer         { return shift->_generic_get('alignment_namer');           }
sub _get_control_alignment_namer { return shift->_generic_get('control_alignment_namer');   }
sub _get_samtools_fasta_index    { return shift->_generic_get('samtools_fasta_index');      }
sub _get_chromosome_lengths_by_species_assembly { 
  return shift->_generic_get('chromosome_lengths_by_species_assembly');
}

sub construct {
  my $self  = shift;
  
  (
    my $signal_alignment_for_peak_calling, 
    my $control_alignment_for_peak_calling
  ) = $self->_construct_bam_to_bed_conversion_plan;
  
  $self->_construct_peak_calling_plan(
    $signal_alignment_for_peak_calling,
    $control_alignment_for_peak_calling
  );
  return;
}

=head2 _construct_bam_to_bed_conversion_plan

  ccat only supports bed format, so if broad peaks have to be called, the 
  bam files must be converted to bed format.

=cut
sub _construct_bam_to_bed_conversion_plan {
  my $self  = shift;

  my $experiment = $self->_get_experiment;
  
  my $peak_calling_strategy 
    = $self->select_peak_calling_strategy(
        $experiment
    );
    
  my $bam_file_to_bed_conversion_plan = [];
  
  my $bam_file_to_bed_conversion_needed 
    = $peak_calling_strategy eq CALL_BROAD_PEAKS;
  
  my $signal_alignment_for_peak_calling  = $self->_get_signal_alignment;
  my $control_alignment_for_peak_calling = $self->_get_control_alignment;
  
  if ($bam_file_to_bed_conversion_needed) {
      my $alignment_namer         = $self->_get_alignment_namer;
      my $control_alignment_namer = $self->_get_control_alignment_namer;
      
      my $signal_alignment_as_bed = {
        input      => $signal_alignment_for_peak_calling,
        type       => ALIGNMENT_ANALYSIS,
        name       => $alignment_namer->base_name_no_duplicates,
        task       => CONVERT_BAM_TO_BED_ANALYSIS,
        is_control => FALSE,
        output     => {
          real => $alignment_namer->bed_file_no_duplicates,
          format => BED_FORMAT
        },
      };
      
      my $control_alignment_as_bed;
      
      if ($experiment->has_control_Experiment) {
        $control_alignment_as_bed = {
            input      => $control_alignment_for_peak_calling,
            type       => ALIGNMENT_ANALYSIS,
            name       => $control_alignment_namer->base_name_no_duplicates,
            task       => CONVERT_BAM_TO_BED_ANALYSIS,
            is_control => TRUE,
            output     => {
            real => $control_alignment_namer->bed_file_no_duplicates,
            format => BED_FORMAT
            },
        };
      }
      
    $bam_file_to_bed_conversion_plan = [
      $signal_alignment_as_bed,
      $control_alignment_as_bed
    ];
    
    $signal_alignment_for_peak_calling  = $signal_alignment_as_bed;
    $control_alignment_for_peak_calling = $control_alignment_as_bed;
  }
  $self->_set_bam_to_bed_conversion_plan($bam_file_to_bed_conversion_plan);
  return (
    $signal_alignment_for_peak_calling, 
    $control_alignment_for_peak_calling
  );
}

sub _construct_peak_calling_plan {
  my $self  = shift;
  
  my $signal_alignment_for_peak_calling  = shift;
  my $control_alignment_for_peak_calling = shift;
  
  my $experiment = $self->_get_experiment;
  
  my $peak_calling_strategy 
    = $self->select_peak_calling_strategy(
        $experiment
    );
  
  my $analysis_logic_name;
  
  my $idr = $self->_get_idr;
  
  my $idr_strategy = $idr->{strategy};
  
  if ($peak_calling_strategy eq CALL_BROAD_PEAKS) {
    $analysis_logic_name = ENSEMBL_BROAD_PEAK_CALLING_ANALYSIS;
  }
  
  # narrow peaks + no idr => peak calling for narrow peaks
  if (
    ($peak_calling_strategy eq CALL_NARROW_PEAKS) && ($idr_strategy eq SKIP_IDR)
  ) {
    $analysis_logic_name = ENSEMBL_NARROW_PEAK_CALLING_ANALYSIS_DEFAULT;
  }
  
  # narrow peaks + idr => permissive peak calling
  if (
    ($peak_calling_strategy eq CALL_NARROW_PEAKS) && ($idr_strategy ne SKIP_IDR)
  ) {
    $analysis_logic_name = ENSEMBL_NARROW_PEAK_CALLING_ANALYSIS_PERMISSIVE;
  }
  
  # tight peaks + no idr => peak calling for tight peaks
  if (
    ($peak_calling_strategy eq CALL_TIGHT_PEAKS) && ($idr_strategy eq SKIP_IDR)
  ) {
    $analysis_logic_name = ENSEMBL_TIGHT_PEAK_CALLING_ANALYSIS_DEFAULT;
  }
  
  # tight peaks + idr => permissive peak calling
  if (
    ($peak_calling_strategy eq CALL_TIGHT_PEAKS) && ($idr_strategy ne SKIP_IDR)
  ) {
    $analysis_logic_name = ENSEMBL_NARROW_PEAK_CALLING_ANALYSIS_PERMISSIVE;
  }

  
  if (! defined $analysis_logic_name) {
    confess("Couldn't assign a logic name for the peak calling strategy $peak_calling_strategy!");
  }
  
  my $feature_type = $experiment->feature_type;
  my $feature_type_name = $feature_type->name;
  
  my $epigenome = $experiment->epigenome;
  my $epigenome_name = $epigenome->name;
  
  my $experimental_group = $experiment->experimental_group;
  my $experimental_group_name = $experimental_group->name;
  
  my $display_label = "$feature_type_name - $epigenome_name enriched sites";
  
  use Bio::EnsEMBL::Funcgen::Utils::GoodUtils qw( create_production_name );
  my $name = create_production_name(
    join " ", 
      $epigenome_name,
      $feature_type_name,
      $analysis_logic_name,
      $experimental_group_name
  );
  
  my $output_dir = $self->_get_peaks_output_dir;
  
  my $peak_calling_plan = {
      input => {
        signal  => $signal_alignment_for_peak_calling,
        control => $control_alignment_for_peak_calling,
      },
      run_idr               => $self->_get_idr,
      peak_calling_strategy => $peak_calling_strategy,
      type                  => PEAK_CALLING_TYPE,
      task                  => CALL_PEAKS_ANALYSIS,
      feature_type          => $feature_type_name,
      epigenome             => $epigenome_name,
      analysis              => $analysis_logic_name,
      output_dir            => $output_dir,
      display_label         => $display_label,
      name                  => $name,
      samtools_fasta_index  => $self->_get_samtools_fasta_index,
      chromosome_lengths_by_species_assembly
        => $self->_get_chromosome_lengths_by_species_assembly,
    }
  ;
  $self->_set_peakcalling_plan($peak_calling_plan);
  return;
}

sub select_peak_calling_strategy {

  my $self = shift;
  my $experiment = shift;
  
  my $feature_type = $experiment->feature_type;

  if ($feature_type->_creates_broad_peaks) {
    return CALL_BROAD_PEAKS;
  }
  if ($feature_type->name eq 'DNase1') {
    return CALL_TIGHT_PEAKS;
  }
  return CALL_NARROW_PEAKS
}

1;
